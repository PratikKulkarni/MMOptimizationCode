/*!
  \file cltest.c
  \brief Test program to exercise the OpenCL API

  This file is a test program to compile and execute kernels
  to exercise and validate the OpenCL API implementation.
*/

/* Include files */
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include<time.h>

//#include "redefine.h"
#define CL_TARGET_OPENCL_VERSION 220
#include "CL/opencl.h"
#include "Gemm_JPI.c"
#define PROGRAM_FILE "mm_dist_blis_opt_final28.c"
#define __NUMCOL__ 7
#define __NUMROW__ 5
#define dabs( x ) ( (x) < 0 ? -(x) : x )

const int M = 128;
const int N = 128;
const int K = 128;
#define MAX_SOURCE_SIZE (0x100000)


float MaxAbsDiff( int m, int n, float *ap, int lda, float *bp, int ldb )
/*
   MaxAbsDiff returns the maximum absolute difference over
   corresponding elements of matrices A and B.
*/
{
  float diff=0.0;
  int  i, j;

  for ( i=0; i<m; i++ )
    for ( j=0; j<n; j++ )
      if ( dabs( ap[i*lda+j] - bp[i*ldb+j] ) > diff )
	  diff = dabs( ap[i*lda+j] - bp[i*ldb+j] );

  return diff;
}


void RandomMatrix( int m, int n, float *ap, int lda )
/*
   RandomMatrix overwrite A with random values.
*/
{
  int  i, j;

	srand(time(NULL));
  for ( i=0; i<m; i++ )
    for ( j=0; j<n; j++ )
      ap[i*lda+j] = (rand() % 10);
}


void checkErr(cl_int err, const char *name){
	if(err != CL_SUCCESS) {
		printf("ERROR: %s\n",name);
    exit(EXIT_FAILURE);
  }
}

void initMatrix(float *a, int n, int m){
  int i = 0;
  int j = 0;
  for(i = 0; i < n; i = i+1){
		for(j = 0; j < m; j = j+1){
			a[i*m+j] = 10 * ((float)rand()/RAND_MAX);
    }
  }
}

void zeroInitMatrix(float *a, int n, int m){
	int i = 0;
	int j = 0;
  for(i = 0; i < n; i = i+1){
		for(j = 0; j < m; j = j+1){
			a[i*m+j] = (float)0;
		}
	}
}

void showMatrix(float *a, int m, int n){
  for (int i = 0; i < m; i++) {
 		for (int j = 0; j < n; j++) {
 			printf("%f ", a[i*n+j]);
 		}
 		printf("\n");
 	}
}

int runmatmul()
{

  float *A;
	float *B;
	float *C;
  float *Atilde;
	float *Btilde;
	float *Cref;
	float *temp1;
	float *temp2;
  float *temp3;
	int *data;
	int szA, szB, szC,val1,val2,iNode[__NUMROW__][__NUMCOL__],jNode[__NUMROW__][__NUMCOL__],x,y;
  int rows = 0;
	int columns,n,m,k,ldA,ldB,ldC;
	float n1,diff,d_one = 1.0;
	int count = 0;
  char ch;
	n=N;
	m=M;
	k=K;
	ldA=K;
	ldB=N;
	ldC=N;
  #define KC 40
  #define MC 130
  #define NC 16
  #define MR 4
  #define NR 4
//  #define MPC 4
  szA  = M*K;
	szB  = K*N;
	szC  = M*N;
	A = (float *)malloc(szA*sizeof(float));
	B = (float *)malloc(szB*sizeof(float));
	C = (float *)malloc(szC*sizeof(float));
	Cref=(float *)malloc(szC*sizeof(float));
	printf("Reached Here 00\n");
	RandomMatrix( m, k, A, ldA );
	//showMatrix(A, m,k);
	RandomMatrix( k, n, B, ldB );
	//showMatrix(B, k,n);
	for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          Cref[i*ldC+j] = 0;
        }
      }
	MyGemm(m,n,k,A,ldA,B,ldB,Cref,ldC);
  int d0,d1;
  d0=__NUMROW__;
  d1=__NUMCOL__;
	if((M/d0+(d0-1))%MC==0)
	val1=(M/d0+(d0-1))/MC;
	else
	val1=(M/d0+(d0-1))/MC+1;
	temp1 = (float *)malloc((M/d0+(d0-1))*K*d0*sizeof(float));
	temp2 = (float *)malloc(K*(N/d1+(d1-1))*d1*sizeof(float));
  temp3 = (float *)malloc((M/d0+(d0-1))*(N/d1+(d1-1))*d0*d1*sizeof(float));
	data = (int *)malloc(val1*11*d0*d1*sizeof(int));
  Atilde = (float *)malloc(MC*KC*val1*d0*d1*sizeof(float));
  Btilde = (float *)malloc(KC*(N/d1+(d1-1))*val1*d0*d1*sizeof(float));
  for(int i=0;i<(M/d0+(d0-1))*d0;i++){
    for(int j=0;j<(N/d1+(d1-1))*d1;j++)
    temp3[i*(N/d1+(d1-1))*d1+j]=0;
  }
  for(int i=0;i<__NUMROW__;i++){
    for(int j=0;j<__NUMCOL__;j++){
      iNode[i][j]=0;
      jNode[i][j]=0;
    }
  }
  for(int i=0;i<M;i++){
    x=(i%d0);

    for(int j=0;j<K;j++){
        temp1[x*(M/d0+(d0-1))*K+iNode[x][0]*K+jNode[x][0]]=A[i*K+j];
        jNode[x][0]++;

    }

    for(int p=0;p<__NUMROW__;p++){
      for(int q=0;q<__NUMCOL__;q++){
        if(p==x)
        iNode[p][q]++;
        jNode[p][q]=0;
      }
    }

  }

  for(int i=0;i<__NUMROW__;i++){
    for(int j=0;j<__NUMCOL__;j++){
      iNode[i][j]=0;
      jNode[i][j]=0;
    }
  }

  /*for(int i=0;i<K;i++){
    for(int j=0;j<(N/d1+(d1-1))*d1;j++){
      temp2[i*(N/d1+(d1-1))*d1+j]=-1;
    }
  }*/

  for(int j=0;j<N;j++){
    y=(j%d1);


    for(int i=0;i<K;i++){

        temp2[y*(N/d1+(d1-1))*K+iNode[0][y]*(N/d1+(d1-1))+jNode[0][y]]=B[i*N+j];
        iNode[0][y]++;

    }
    for(int p=0;p<__NUMROW__;p++){
      for(int q=0;q<__NUMCOL__;q++){
        iNode[p][q]=0;
        if(q==y)
        jNode[p][q]++;
      }
    }
  }



  for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          C[i*ldC+j] = 0;
        }
      }



  cl_int error = 0;
  // INITIALIZATION

  cl_context contextId = NULL;
  cl_command_queue cmdqId = NULL;
  cl_device_id deviceId = NULL;

  cl_device_id *devices;
  cl_uint num_devices;
  char name_data[48];
  cl_platform_id platform_id[3];
  cl_uint ret_num_platforms;
  cl_int ret;
  cl_int err;

  /* Identify a platform */
  /* Get Platform and Device Info */
  ret = clGetPlatformIDs(3, platform_id, &ret_num_platforms);
  checkErr(ret,"platform_id");
  printf("No. of platforms detected: %d\n",ret_num_platforms);

  cl_platform_id pocl_platform_id;
  char* pocl_platform_name = "Portable Computing Language";
  for(int i=0; i<ret_num_platforms; i++){
    size_t namelen;
    ret = clGetPlatformInfo(platform_id[i],CL_PLATFORM_NAME,1024,NULL,&namelen);
    char* pname = calloc(namelen, sizeof(char));
    ret = clGetPlatformInfo(platform_id[i],CL_PLATFORM_NAME,namelen,pname,NULL);
    checkErr(ret,"platform_info");
    printf("platform name: %s\n",pname);
    if(strcmp(pname, pocl_platform_name) == 0){
      printf("FOUND\n");
      pocl_platform_id = platform_id[i];
    }
  }


  err = clGetDeviceIDs(pocl_platform_id, CL_DEVICE_TYPE_ACCELERATOR,2, NULL, &num_devices);
	printf("Number of devices : %d\n",num_devices);

  devices = (cl_device_id*) malloc(sizeof(cl_device_id) * num_devices);
  clGetDeviceIDs(pocl_platform_id, CL_DEVICE_TYPE_ACCELERATOR, num_devices, devices, NULL);

  for(cl_uint i=0; i<num_devices; i++) {
    err = clGetDeviceInfo(devices[i], CL_DEVICE_NAME, sizeof(name_data), name_data, NULL);
    if(err < 0) {
      perror("Couldn't read name data");
      exit(1);
    }
    if(strcmp(name_data,"rsim-hyperop")==0){
      printf("name_data : %s\n",name_data);
      deviceId = devices[i];
      break;
    }

  }


  contextId = clCreateContext(NULL, 1, &deviceId, NULL, NULL, &error);
	if (error != CL_SUCCESS){
		if (error == CL_INVALID_PLATFORM)
		printf("CL_INVALID_PLATFORM\n");
		else if (error == CL_INVALID_PROPERTY)
		printf("CL_INVALID_PROPERTY\n");
		else if (error == CL_INVALID_VALUE)
		printf("CL_INVALID_VALUE\n");
		else if (error == CL_DEVICE_NOT_AVAILABLE)
		printf("CL_DEVICE_NOT_AVAILABLE\n");
		else if (error == CL_OUT_OF_RESOURCES)
		printf("CL_OUT_OF_RESOURCES\n");
		else if (error == CL_OUT_OF_HOST_MEMORY)
		printf("CL_OUT_OF_HOST_MEMORY\n");
		else if (error == CL_INVALID_DEVICE)
		printf("CL_INVALID_DEVICE\n");
		else if (error == CL_INVALID_OPERATION)
		printf("CL_INVALID_OPERATION\n");
		else if (error == CL_INVALID_OPERATION)
		printf("CL_INVALID_OPERATION\n");
	}
  checkErr(error, "contextId");

	printf("Checkpoint0\n");

  cmdqId = clCreateCommandQueueWithProperties( contextId, deviceId, NULL , &error);
  checkErr(error, "commandQueueId");

	printf("Checkpoint1\n");

  //CREATE KERNEL
  FILE *program_handle;
  char *program_buffer;
  size_t program_size;
  cl_program programId;

  program_handle = fopen(PROGRAM_FILE, "r");
  if(program_handle == NULL) {
     perror("Couldn't find the program file");
     exit(1);
  }
  fseek(program_handle, 0, SEEK_END);
  program_size = ftell(program_handle);
  rewind(program_handle);
  program_buffer = (char*)malloc(program_size + 1);
  program_buffer[program_size] = '\0';
  fread(program_buffer, sizeof(char), program_size, program_handle);
  fclose(program_handle);

  /* Create program from file */
  programId = clCreateProgramWithSource(contextId, 1,
     (const char**)&program_buffer, &program_size, &error);
  if(error < 0) {
     perror("Couldn't create the program");
     exit(1);
  }
  free(program_buffer);

  printf("Reached here 3\n");

  const char *compile_options = {"-ffp-contract=fast"};
  error = clBuildProgram(programId, 1, &deviceId, compile_options, NULL, NULL);

  if (error != CL_SUCCESS){
    char buildLog[16384];
    clGetProgramBuildInfo(programId, deviceId, CL_PROGRAM_BUILD_LOG , sizeof(buildLog) , buildLog , NULL);
    printf("Buffer : %s\n", buildLog);
  	return 1;
  }

  checkErr(error, "buildProgram");

  printf("Reached here 4\n");

  cl_kernel kernelId;
  kernelId = clCreateKernel(programId, "matmulstart", &error);
  checkErr(error, "kernelId");

  //CREATE KERNEL I/O BUFFERS

  cl_mem dev_temp1,dev_temp2,dev_temp3,dev_data,dev_Atilde,dev_Btilde;


	dev_temp1 = clCreateBuffer(contextId, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, (M/d0+(d0-1))*K*d0*sizeof(float), temp1, &error);
	checkErr(error, "CreateBuffer dev_temp1");

	dev_temp2 = clCreateBuffer(contextId, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, K*(N/d1+(d1-1))*d1*sizeof(float), temp2, &error);
	checkErr(error, "CreateBuffer dev_temp2");

  dev_temp3 = clCreateBuffer(contextId, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, (M/d0+(d0-1))*(N/d1+(d1-1))*d0*d1*sizeof(float), temp3, &error);
	checkErr(error, "CreateBuffer dev_temp3");

  dev_Atilde = clCreateBuffer(contextId, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, MC*KC*val1*d0*d1*sizeof(float), Atilde, &error);
	checkErr(error, "CreateBuffer dev_Atilde");

	dev_Btilde = clCreateBuffer(contextId, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, KC*(N/d1+(d1-1))*val1*d0*d1*sizeof(float), Btilde, &error);
	checkErr(error, "CreateBuffer dev_Btilde");

	dev_data = clCreateBuffer(contextId, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, val1*11*d0*d1*sizeof(int), data, &error);
	checkErr(error, "CreateBuffer dev_data");

  //SET KERNEL ARGUMENTS

	error = clSetKernelArg(kernelId, 0, sizeof(cl_mem), &dev_temp1); checkErr(error, "SetKernelArg dev_temp1");
	error = clSetKernelArg(kernelId, 1, sizeof(cl_mem), &dev_temp2); checkErr(error, "SetKernelArg dev_temp2");
  error = clSetKernelArg(kernelId, 2, sizeof(cl_mem), &dev_temp3); checkErr(error, "SetKernelArg dev_temp3");
  error = clSetKernelArg(kernelId, 3, sizeof(cl_mem), &dev_Atilde); checkErr(error, "SetKernelArg dev_Atilde");
  error = clSetKernelArg(kernelId, 4, sizeof(cl_mem), &dev_Btilde); checkErr(error, "SetKernelArg dev_Btilde");
	error = clSetKernelArg(kernelId, 5, sizeof(cl_mem), &dev_data); checkErr(error, "SetKernelArg dev_data");
  //int nDim = ARRAY_SIZE;
  error = clSetKernelArg(kernelId, 6, sizeof(int), &m); checkErr(error, "SetKernelArg &m");
  error = clSetKernelArg(kernelId, 7, sizeof(int), &n); checkErr(error, "SetKernelArg &n");
  error = clSetKernelArg(kernelId, 8, sizeof(int), &k); checkErr(error, "SetKernelArg &k");

  //KERNEL EXECUTION LOOP

	printf("Reached here 5\n");
  //SEND INPUTS


	error = clEnqueueWriteBuffer(cmdqId, dev_temp1, 1, 0, (M/d0+(d0-1))*K*d0*sizeof(float), temp1, 0, NULL, NULL);
  checkErr(error, "WriteBuffer dev_temp1");


	error = clEnqueueWriteBuffer(cmdqId, dev_temp2, 1, 0, K*(N/d1+(d1-1))*d1*sizeof(float), temp2, 0, NULL, NULL);
  checkErr(error, "WriteBuffer dev_temp2");


  error = clEnqueueWriteBuffer(cmdqId, dev_temp3, 1, 0, (M/d0+(d0-1))*(N/d1+(d1-1))*d0*d1*sizeof(float), temp3, 0, NULL, NULL);
  checkErr(error, "WriteBuffer dev_temp3");

  error = clEnqueueWriteBuffer(cmdqId, dev_Atilde, 1, 0,MC*KC*val1*d0*d1*sizeof(float), Atilde, 0, NULL, NULL);
  checkErr(error, "WriteBuffer dev_temp3");

  error = clEnqueueWriteBuffer(cmdqId, dev_Btilde, 1, 0, KC*(N/d1+(d1-1))*val1*d0*d1*sizeof(float), Btilde, 0, NULL, NULL);
  checkErr(error, "WriteBuffer dev_temp3");


	error = clEnqueueWriteBuffer(cmdqId, dev_data, 1, 0, val1*11*d0*d1*sizeof(int), data, 0, NULL, NULL);
  checkErr(error, "WriteBuffer dev_data");

	printf("Reached here 6\n");

  size_t global[1] = {1};
  size_t local[1] = {1};

	printf("Reached here 7\n");

  //LAUNCH KERNEL
  error = clEnqueueNDRangeKernel(cmdqId, kernelId, 1, NULL, global, local, 0, NULL, NULL);
	if (error != CL_SUCCESS){
		if (error == CL_INVALID_PROGRAM_EXECUTABLE)
		printf("CL_INVALID_PROGRAM_EXECUTABLE\n");
		else if (error == CL_INVALID_COMMAND_QUEUE)
		printf("CL_INVALID_COMMAND_QUEUE\n");
		else if (error == CL_INVALID_KERNEL)
		printf("CL_INVALID_KERNEL\n");
		else if (error == CL_INVALID_CONTEXT)
		printf("CL_INVALID_CONTEXT\n");
		else if (error == CL_INVALID_KERNEL_ARGS)
		printf("CL_INVALID_KERNEL_ARGS\n");
		else if (error == CL_INVALID_WORK_DIMENSION)
		printf("CL_INVALID_WORK_DIMENSION\n");
		else if (error == CL_INVALID_WORK_GROUP_SIZE)
		printf("CL_INVALID_WORK_GROUP_SIZE\n");
		else if (error == CL_INVALID_GLOBAL_OFFSET)
		printf("CL_INVALID_GLOBAL_OFFSET\n");
		else if (error == CL_OUT_OF_RESOURCES)
		printf("CL_OUT_OF_RESOURCES\n");
		else if (error == CL_MEM_OBJECT_ALLOCATION_FAILURE)
		printf("CL_MEM_OBJECT_ALLOCATION_FAILURE\n");
		else if (error == CL_INVALID_EVENT_WAIT_LIST)
		printf("CL_INVALID_EVENT_WAIT_LIST\n");
		else if (error == CL_OUT_OF_HOST_MEMORY)
		printf("CL_OUT_OF_HOST_MEMORY\n");
	}
  checkErr(error, "KernelLaunch");

	printf("Reached here 8\n");



  //READ OUPUTS
  error = clEnqueueReadBuffer(cmdqId, dev_temp3, 1, 0, (M/d0+(d0-1))*(N/d1+(d1-1))*d0*d1*sizeof(float), temp3, 0, NULL, NULL);
  checkErr(error, "ReadBuffer");


  for(int i=0;i<__NUMROW__;i++){
    for(int j=0;j<__NUMCOL__;j++){
      x=0;
      for(int p=i;p<m;p+=d0){
        y=0;
        for(int q=j;q<n;q+=d1){
          C[p*n+q]=temp3[i*(M/d0+(d0-1))*(N/d1+(d1-1))*d1+j*(M/d0+(d0-1))*(N/d1+(d1-1))+x*(N/d1+(d1-1))+y];
          y++;
        }
        x++;
      }
    }
  }
  //printf("Matrix C :\n");
  //showMatrix(C, m, n);
  //printf("Matrix Cref :\n");
  //showMatrix(Cref, m, n);
	diff = MaxAbsDiff( m, n, C, ldC, Cref, ldC );
	printf( "%8.4le\n", diff  );


  printf("ALL OK WE ARE DONE\n");


  //RELEASE DEVICE MMEORY OBJECTS

	error = clReleaseMemObject(dev_temp1); checkErr(error, "Release dev_temp1");
  error = clReleaseMemObject(dev_temp2); checkErr(error, "Release dev_temp2");
  error = clReleaseMemObject(dev_temp3); checkErr(error, "Release dev_temp3");
  error = clReleaseMemObject(dev_Atilde); checkErr(error, "Release dev_Atilde");
  error = clReleaseMemObject(dev_Btilde); checkErr(error, "Release dev_Btilde");
	error = clReleaseMemObject(dev_data); checkErr(error, "Release dev_data");


  //---------
  error = clReleaseKernel(kernelId); checkErr(error, "ReleaseKernel");
  error = clReleaseProgram(programId); checkErr(error, "ReleaseProgram");

  error = clReleaseCommandQueue(cmdqId); checkErr(error, "ReleaseCommandQueue");
  error = clReleaseContext(contextId); checkErr(error, "ReleaseContextId");
  printf("Reached here 2\n");
  return 0;
}

int main(int argc, char* argv[])
{

/*  int n = 32;
  if(argc > 1){
    n = atoi(argv[1]);
    if(n%8 != 0){
      printf("INFO:Matrix size must be a multiple of 8, changing it to 32\n");
      n = 32;
    }
  }*/



  printf("#####################################################\n");
  printf("##########  TESTING MATRIX MULTIPLICATION  ##########\n");
  printf("#####################################################\n");

  runmatmul();

}
