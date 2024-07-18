import subprocess

def load_docker_image(image_tarball_path):
    try:
        # Construct the docker load command
        docker_load_command = ["docker", "load", "-i", image_tarball_path]
        
        # Execute the docker load command
        subprocess.run(docker_load_command, check=True)
        
        print("Docker image loaded successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error: Failed to load Docker image. {e}")

def run_docker_command(command):
    try:
        # Execute the docker command
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: Failed to run Docker command. {e}")

def main():
    # Prompt the user for the path to the Docker image tarball file
    image_tarball_path = r"C:\Users\ferra\Documents\codici\CH4\CTMF\ctmf\ctmf_v1.0.tar"
    
    # Load the Docker image from the tarball
    load_docker_image(image_tarball_path)

    # Specify the Docker command with all necessary inputs
    docker_command = """docker run -v D:/CLEAR_UP/CH4_detection:/app/data ctmf \
                     "/app/data/img_PRISMA/USA/Encinal_wells_USA/20200426/PRS_L1_STD_OFFL_20200426172708_20200426172712_0001.he5" \
                     "/app/data/img_PRISMA/USA/Encinal_wells_USA/20200426/PRS_L2C_STD_20200426172708_20200426172712_0001.he5" \
                     "/app/data/Matched_filter_approach/codici/CTMF/DEM_1Km/srtm30plus_v11_land.nc" \
                     "/app/data/Matched_filter_approach/codici/CTMF/LUTs/dataset_ch4_full.hdf5" \
                     "/app/data/img_PRISMA/USA/Encinal_wells_USA/20200426/MF_out" \
                     --clusters 30"""
    
    # Run the Docker command
    run_docker_command(docker_command)

if __name__ == "__main__":
    main()
