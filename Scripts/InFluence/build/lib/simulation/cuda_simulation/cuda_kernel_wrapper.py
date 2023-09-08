import os

def get_script_directory():
    # Get the directory where the currently executing script is located
    script_path = os.path.abspath(__file__)
    script_dir = os.path.dirname(script_path)
    return script_dir

def load_cuda_kernel_code(file_path):
    with open(file_path, "r") as f:
        return f.read()

def get_combined_cuda_kernel_code():
    # Get the directory of the currently executing script
    script_dir = get_script_directory()

    # Define the names of CUDA kernel files relative to the script directory
    kernel_files = [
        "gpu_helpers.cu",
        "MCLoop.cu"
    ]

    combined_code = ""
    for file_name in kernel_files:
        file_path = os.path.join(script_dir, file_name)
        combined_code += load_cuda_kernel_code(file_path) + "\n"

    return combined_code


