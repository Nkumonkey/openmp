#include<iostream>
#include<arm_neon.h>
#include <omp.h>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#define NUM_THREADS 7
using namespace std;
int n;
 float **A;
void f_omp_static_neon ()
{
float32x4_t va = vmovq_n_f32(0) ;
float32x4_t vx = vmovq_n_f32(0) ;
float32x4_t vaij = vmovq_n_f32(0) ;
float32x4_t vaik = vmovq_n_f32(0) ;
float32x4_t vakj = vmovq_n_f32(0) ;
#pragma omp p a r a l l e l num_threads (NUM_THREADS) , private (va , vx , vaij , vaik, vakj )
for (int k = 0; k < n ; k++){
#pragma omp s i n g l e //串 行 部 分
{
float32x4_t vt=vmovq_n_f32(A[ k ] [ k ] ) ; int j ;
for ( j = k + 1; j < n ; j++){
    va=vld1q_f32(&(A[ k ] [ j ] ) ) ;
va= vdivq_f32 (va , vt ) ;
vst1q_f32(&(A[ k ] [ j ] ) , va ) ;
}
for ( ; j<n ; j++)
A[ k ] [ j ]=A[ k ] [ j ] * 1.0 / A[ k ] [ k ] ;
A[ k ] [ k ] = 1.0 ;
}
#pragma omp for schedule ( static ) //并 行 部 分
for (int i = k + 1; i < n ; i++){
vaik=vmovq_n_f32(A[ i ] [ k ] ) ; int j ;
for ( j = k + 1; j+4 <= n ; j+=4){
vakj=vld1q_f32(&(A[ k ] [ j ] ) ) ;
vaij=vld1q_f32(&(A[ i ] [ j ] ) ) ;
vx=vmulq_f32 ( vakj , vaik ) ;
vaij=vsubq_f32 ( vaij , vx) ;
vst1q_f32(&A[ i ] [ j ] , vaij ) ;
}
for ( ; j<n ; j++)
A[ i ] [ j ] = A[ i ] [ j ] - A[ i ] [ k ] * A[ k ] [ j ] ;
A[ i ] [ k ] = 0;
}
}
}

 int main(){
     cin>>n;
      srand(time(NULL));
      for(int i=0;i<n;++i)
        A=new float* [n];
     for(int i=0;i<n;++i)
        A[i]=new float[n]();
     for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i][j] = rand(); // 生成一个 0 到 9 的随机数
        }
    }
    struct timeval start, end;
    double interval;
    gettimeofday(&start, NULL);
    f_omp_static_neon ();
     gettimeofday(&end, NULL);
      interval = 1000000*(end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec);
 for(int i=0;i<n;++i){
 for(int j=0;j<n;++j)
 cout<<A[i][j]<<" ";
 cout<<endl;
 }
 cout<<interval<<endl;
 return 0;
 }
