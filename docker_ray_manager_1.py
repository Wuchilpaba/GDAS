import docker,os,time
import resource_path

client = docker.from_env()
current_directory = os.getcwd()
worker_ray_path = resource_path.resource_path('resources')
def start_ray_container():
    print("Starting Ray container...")
    time.sleep(3)
    container = client.containers.run(
        image="bitnami/ray:latest",
        detach=True,
        tty=True,
        auto_remove=True,
        volumes={
            worker_ray_path: {
                "bind": "/app",
                "mode": "rw"
            }
        },
        command=f"/app/task_1.py",
        name="bootstrap_iteration_By"
    )
    print(f"Started container:{container.id}")
    return container

def stop_ray_container(container):
    task_complete_flag = "/app/task1_complete.txt"
    while True:
        if container.exec_run(f"test -f {task_complete_flag}").exit_code == 0:
            print(f"Task in container {container.id} completed successfully.")
            break
        else:
            print(f"Waiting for task in container {container.id} to complete...")
            time.sleep(2)
    print("Stopping Ray container...")
    container.stop()
    while not os.path.join(worker_ray_path,'task1_complete.txt'):
        time.sleep(2)
        print(f"Waiting for task in container {container.id} to stop...")