import pickle
import ray,numpy as np, os

@ray.remote
def bootstrap_iteration_By(data1, data2):
    # Perform one bootstrap sample iteration
    sample1 = np.random.choice(data1, size=len(data1), replace=True)
    sample2 = np.random.choice(data2, size=len(data2), replace=True)
    return np.mean(sample1) - np.mean(sample2)

def main():
    current_directory = os.getcwd()
    worker_ray_path = os.path.join(current_directory)
    ray.init(ignore_reinit_error=True)
    with open(rf'{worker_ray_path}/data_cache1.pickle', 'rb') as f:
        data = pickle.load(f)
    group1_data_scaled = data[0]
    group2_data_scaled = data[1]
    futures = [bootstrap_iteration_By.remote(group1_data_scaled , group2_data_scaled) for _ in range(3000)]
    result = ray.get(futures)
    print(result)
    with open(rf'{worker_ray_path}/result_cache1.pickle', 'wb') as f:
        pickle.dump(result, f)
    with open("/app/task1_complete.txt", "w") as f:
        f.write("Task1 completed successfully")
    print('Task1 Completed')

if __name__ == "__main__":
    main()
