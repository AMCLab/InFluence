# cuda_kernels_wrapper.py

def load_cuda_kernel_code(file_path):
    with open(file_path, "r") as f:
        return f.read()

def get_combined_cuda_kernel_code():
    # Load individual CUDA kernel files and concatenate them into a single string
    kernel_files = [
        "/home/uclworkstation1/InFluence/Scripts/InFluence/CUDA_FILE/MCLoop.cu"
        
        
        
    ]

    combined_code = ""
    for file_path in kernel_files:
        combined_code += load_cuda_kernel_code(file_path) + "\n"

    return combined_code
