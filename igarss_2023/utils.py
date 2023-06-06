import shutil
import subprocess
from datetime import datetime, timedelta
from pathlib import Path

import asf_search as asf
import hyp3_sdk as sdk
import pandas as pd
import requests


SCENE_LIST_URL = 'https://zenodo.org/record/7151431/files/InSAR_scenes.txt'


def download_to_file(url, filename):
    if Path(filename).exists():
        return

    with open(filename, 'wb') as f:
        response = requests.get(url)
        f.write(response.content)


def get_pair_list(lines, direction, orbit):
    pairs = [pair.split('_') for pair in lines]
    dates = sorted(list({datetime.strptime(item, '%Y%m%d') for sublist in pairs for item in sublist}))
    direction = asf.ASCENDING if direction == 'ascending' else asf.DESCENDING
    date_pre = (min(dates) - timedelta(days=1)).strftime('%Y-%m-%d')
    date_post = (max(dates) + timedelta(days=1)).strftime('%Y-%m-%d')
    search_results = asf.geo_search(
        platform=asf.SENTINEL1,
        intersectsWith='POINT(-135.7525 57.0524)',
        start=date_pre,
        end=date_post,
        processingLevel=asf.SLC,
        beamMode=asf.IW,
        flightDirection=direction,
        relativeOrbit=orbit,
    )

    names = [item.properties['sceneName'] for item in search_results]
    rows = []
    for date in dates:
        scene = [name for name in names if name[17:25] == date.strftime('%Y%m%d')]
        if len(scene) != 1:
            raise ValueError(f'Found {len(scene)} scenes, should find exactly 1.')
        rows.append([scene[0], date])

    scene_df = pd.DataFrame(rows, columns=['scene', 'date'])

    scene_pairs = []
    for pair in pairs:
        date1 = datetime.strptime(pair[0], '%Y%m%d')
        date2 = datetime.strptime(pair[1], '%Y%m%d')
        scene1 = scene_df.loc[scene_df['date'] == date1, 'scene'].values[0]
        scene2 = scene_df.loc[scene_df['date'] == date2, 'scene'].values[0]
        scene_pairs.append([date1, date2, scene1, scene2])

    scene_df = pd.DataFrame(scene_pairs, columns=['date1', 'date2', 'scene1', 'scene2'])

    start_date = datetime(2018, 1, 1)
    stop_date = datetime(2020, 1, 1)
    smaller_scenes = scene_df.loc[(scene_df['date1'] > start_date) & (scene_df['date2'] < stop_date)]
    smaller_scenes.to_csv(f'{direction.lower()}_{orbit}_pairs.csv', index=False)


def get_edgecumbe_insar_pairs():
    download_to_file(SCENE_LIST_URL, 'InSAR_scenes.txt')
    input_file = Path('InSAR_scenes.txt').read_text().split('\n')
    stacks = {
        'descending_174': [810, 1123],
        'ascending_79': [433, 708],
        'ascending_50': [92, 343],
    }

    for stack, (start, stop) in stacks.items():
        direction, orbit = stack.split('_')
        get_pair_list(input_file[start:stop], direction, int(orbit))


def submit_stack(pair_file, project_name, water_mask=True):
    pair_list = pd.read_csv(pair_file, parse_dates=[0, 1]).sort_values('date1').reset_index(drop=True)
    sbas_pairs = list(zip(list(pair_list['scene1']), list(pair_list['scene2'])))
    hyp3 = sdk.HyP3()
    jobs = sdk.Batch()
    for reference, secondary in sbas_pairs:
        jobs += hyp3.submit_insar_job(
            reference,
            secondary,
            name=project_name,
            include_dem=True,
            include_look_vectors=True,
            apply_water_mask=water_mask,
        )
    return project_name


def download_stack(project_name, download_dir):
    hyp3 = sdk.HyP3()
    jobs = hyp3.find_jobs(name=project_name)
    statuses = [job.status_code == 'SUCCEEDED' for job in jobs]
    if sum(statuses) != len(statuses):
        print('Jobs are not ready, watching the progress')
        hyp3.watch(jobs)
    insar_products = jobs.download_files(download_dir)
    insar_products = [sdk.util.extract_zipped_product(ii) for ii in insar_products]


def prep_stack_for_mintpy(pairs_file, data_dir, mintpy_dir):
    roi = [6315473, 6340109, 446628, 466268]
    reference_point = [6330696, 456350]
    mintpy_config = mintpy_dir / 'mintpy_config.cfg'
    _ = mintpy_config.write_text(
        f"""
    ##---------processor:
    mintpy.load.processor        = hyp3
    mintpy.plot                  = no

    ##---------interferogram datasets:
    mintpy.load.unwFile          = ../{data_dir}/*/*_unw_phase_clipped.tif
    mintpy.load.corFile          = ../{data_dir}/*/*_corr_clipped.tif

    ##---------geometry datasets:
    mintpy.load.demFile          = ../{data_dir}/*/*_dem_clipped.tif
    mintpy.load.incAngleFile     = ../{data_dir}/*/*_lv_theta_clipped.tif
    mintpy.load.azAngleFile      = ../{data_dir}/*/*_lv_phi_clipped.tif
    mintpy.load.waterMaskFile    = ../{data_dir}/*/*_water_mask_clipped.tif

    ##--------dataset geographic subset
    mintpy.subset.lalo           = [{roi[0]}:{roi[1]}, {roi[2]}:{roi[3]}]
    mintpy.reference.lalo        = {reference_point}

    ##--------network selections
    mintpy.network.coherenceBased  = yes
    mintpy.network.minCoherence    = 0.7

    ##--------unwrapping error correction
    mintpy.unwrapError.method      = no

    ##--------troposperic phase correction
    mintpy.troposphericDelay.method = pyaps # pyaps eventually
    mintpy.troposphericDelay.weatherModel = ERA5
    mintpy.troposphericDelay.weatherDir   = ../weather

    ##--------topographic phase correction
    mintpy.topographicResidual     = yes
    """
    )
    shutil.copy(pairs_file, mintpy_dir / pairs_file)
    cmd = f'smallbaselineApp.py --dir {mintpy_dir} {mintpy_config} --dostep load_data'
    subprocess.run(cmd.split(' '), shell=True, check=True)


def main():
    get_edgecumbe_insar_pairs()
    stacks = [
        ('descending_174_pairs.csv', 'descending_174_nomask'),
        ('descending_174_pairs.csv', 'descending_174'),
        ('ascending_79_pairs.csv', 'ascending_79'),
        ('ascending_50_pairs.csv', 'ascending_50'),
    ]

    project_names = []
    for pair_csv, base_name in stacks:
        if 'nomask' in base_name
            use_mask = True
        else:
            use_mask = False

        project_name = base_name + f'_{datetime.now().strftime("%Y%m%dT%H:%M")}'
        print(f'Project Name is: {project_name}')

        submit_stack(pair_csv, project_name, water_mask=use_mask)
        project_names.append((base_name, project_name))

    pd.DataFrame(project_names).to_csv('projects.txt', sep='\t', index=False, header=False)

    project_names = [['descending_174_pairs.csv', 'descending_174', 'descending_174_20230605T14:01']]
    for pair_csv, base_name, project_name in project_names:
        data_dir = Path('.') / f'{base_name}_data'
        data_dir.mkdir(exist_ok=True)
        download_stack(project_name, data_dir)

        mintpy_dir = Path('.') / f'{base_name}_data'
        mintpy_dir.mkdir(exist_ok=True)
        prep_stack_for_mintpy(pair_csv, data_dir, mintpy_dir)
        shutil.make_archive(f'{base_name}.zip', 'zip', mintpy_dir)


if __name__ == '__main__':
    main()
