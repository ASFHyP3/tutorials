import asf_search as asf
import requests
from pathlib import Path
from datetime import datetime, timedelta
import pandas as pd

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
    # scene_df.to_csv(f'{direction}_{orbit}_scenes.csv', index=False)

    scene_pairs = []
    for pair in pairs:
        date1 = datetime.strptime(pair[0], '%Y%m%d')
        date2 = datetime.strptime(pair[1], '%Y%m%d')
        scene1 = scene_df.loc[scene_df['date'] == date1, 'scene'].values[0]
        scene2 = scene_df.loc[scene_df['date'] == date2, 'scene'].values[0]
        scene_pairs.append([date1, date2, scene1, scene2])

    scene_df = pd.DataFrame(scene_pairs, columns=['date1', 'date2', 'scene1', 'scene2'])
    # scene_df.to_csv(f'{direction}_{orbit}_pairs.csv', index=False)

    start_date = datetime(2018, 1, 1)
    stop_date = datetime(2020, 1, 1)
    smaller_scenes = scene_df.loc[(scene_df['date1'] > start_date) & (scene_df['date2'] < stop_date)]
    smaller_scenes.to_csv(f'{direction}_{orbit}_pairs_smaller.csv', index=False)


if __name__ == '__main__':
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
