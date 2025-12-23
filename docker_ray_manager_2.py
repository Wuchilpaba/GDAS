import docker,os,time
from resource_path import resource_path
client = docker.from_env()
current_directory = os.getcwd()
worker_ray_path = resource_path('resources')
def start_ray_container():
    print("Starting Ray container...")
    time.sleep(3)
    container = client.containers.run(
        image="ray-bayesian:latest",
        detach=True,
        tty=True,
        auto_remove=True,
        volumes={
            worker_ray_path: {
                "bind": "/app",
                "mode": "rw"
            }
        },
        command=f"/app/task_2.py",
        name="bootstrap_iteration"
    )
    print(f"Started container:{container.id}")
    return container

def stop_ray_container(container):
    client = docker.from_env()
    try:
        if client.containers.get(container.id):
            while not os.path.join(worker_ray_path, 'task2_complete.txt'):
                time.sleep(2)
                print(f"Waiting for task in container {container.id} to stop...")
            print("Stopping Ray container...")
            container.stop()
        else:
            print(f"Ray container has been stopped")
    except Exception as e:
        print(e)