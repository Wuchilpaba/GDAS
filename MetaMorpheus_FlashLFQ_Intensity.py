import docker
from docker.errors import DockerException
import os,pickle
import csv,pathlib
import pandas as pd
import logging
import traceback
from resource_path import resource_path
class MetaMorpheus_FlashLFQ_Intensity:
    @staticmethod
    def read_sequences_from_csv(file_path):
        data = []
        try:
            with open(file_path, 'r') as file:
                reader = csv.reader(file, delimiter='\t')
                header = list(filter(None,next(reader)))
                for row in reader:
                    if len(row) == len(header):
                        seq_dict = dict(zip(header, row))
                        data.append(seq_dict)
        except Exception as e:
            logging.error(f"Failed to read sequences from CSV: {e}")
            print(f"Failed to read sequences from CSV: {e}")
            traceback.print_exc()
        df = pd.DataFrame(data)
        filtered_df = df.dropna(subset=['Plausible GlycanComposition'])
        return filtered_df
    @staticmethod
    def read_sequences_from_tsv(file_path):
        try:
            df = pd.read_csv(file_path, sep='\t')
            condidion1 = abs(df['Peak intensity'] >= 1e-9)
            # condition2 = df['Full Sequence'].str.contains('O-Glycosylation')
            filtered_df = df[condidion1]
            return filtered_df
        except Exception as e:
            logging.error(f"Failed to read sequences from TSV: {e}")
            print(f"Failed to read sequences from TSV: {e}")
            traceback.print_exc()

    @classmethod
    def run(cls, PSM_DataFrame, Peak_DataFrame):
        worker_ray_path = pathlib.Path(resource_path('resources')) / "task_Meta"
        data = [PSM_DataFrame, Peak_DataFrame]
        with open(f'{worker_ray_path}\\data_cache_Meta_Intensity.pkl', 'wb') as f:
            pickle.dump(data, f)
        client = docker.from_env()
        print("Starting Ray container...")
        try:
            container = client.containers.run(
                image="bitnami/ray-with-modin:latest",
                detach=True,
                tty=True,
                auto_remove=True,
                volumes={
                    worker_ray_path: {
                        "bind": "/app",
                        "mode": "rw"
                    }
                },
                command=f"/app/task_Meta.py",
                name="MetaMorpheus_FlashLFQ_Intensity"
            )
            return container
        except DockerException as e:
            print(f"Failed to start Ray container: {e}")

