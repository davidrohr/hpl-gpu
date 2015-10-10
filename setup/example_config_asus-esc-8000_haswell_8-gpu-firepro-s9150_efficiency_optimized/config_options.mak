#Select BLAS Backend to use. Available options: MKL, ACML, GOTOBLAS
BLAS_BACKEND				= MKL

#Select which GPU backends are compiled. (The CPU backend is always compiled)
INCLUDE_OPENCL				= 1
INCLUDE_CAL				= 0
INCLUDE_CUDA				= 0
INCLUDE_CUBLAS				= 0

CONFIG_LTO				= 1

#Mark CONFIGURED = 1 to enable compilation
CONFIGURED				= 1
