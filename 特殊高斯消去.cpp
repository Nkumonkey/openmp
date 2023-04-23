#include <arm_neon.h>
#include <stdlib.h>
#include <mutex>
#include <stdio.h>
#include <iostream>
#include <sys/time.h>
#include <time.h>
using namespace std;
#define NUM_THREADS 7
//线程数据结构定义
typedef struct {
    int t_id; // 线程 id
    int** matr;
    int n;
    int m;
    int l;
}threadParam_t;
pthread_rwlock_t rwlock; //读写锁
int shuru(int* a, int n, int m)
{
    for (int j = 0; j < m; j++) a[j] = 0;
    int q = (rand() % n) / 2;
    while (q == 0) q = (rand() % n) / 2;
    int p = n;
    int temp;
    for (int i = 0; i < q; i++)
    {
        if (p != 0)
        {
            p = rand() % p;
            if (i == 0)
                temp = p;
        }
        int bpos = p / 32, spos = p % 32;
        bpos = m - 1 - bpos;
        int ans = 1;
        ans = ans << spos;
        a[bpos] |= ans;
    }
    return temp;
}
void* Static_threadFunc(void* param)
{
    threadParam_t* p = (threadParam_t*)param;
    int n = p->n;
    int m = p->m;
    int l = p->l;
    int** matr = p->matr;
    int t_id = p->t_id;
    int* a = (int*)aligned_alloc(32, n * sizeof(int));
    for (int k = t_id;k < l;k += NUM_THREADS)
    {
        int temp = shuru(a, n, m);
        for (int i = 0; i < n; i++)
        {
            int bpos = i / 32, spos = i % 32;
            if (a[bpos] >> (31 - spos) != 0)
            {
                if (matr[n - 1 - i][bpos] >> (31 - spos) == 0)
                {
                    pthread_rwlock_wrlock(&rwlock);      //写者加写锁
                    for (int j = 0; j < m; j++)
                        matr[n - 1 - i][j] = a[j];
                    pthread_rwlock_unlock(&rwlock);      //释放写锁
                    break;
                }
                else
                {
                    pthread_rwlock_rdlock(&rwlock);    //读者加读锁
                    for (int j = 0; j < m; j++)
                        a[j] ^= matr[n - 1 - i][j];
                    pthread_rwlock_unlock(&rwlock);    //读者释放读锁
                }
            }
        }
    }
    pthread_exit(NULL);
}
void Static_Pthread(int n, int m, int l, int** matr)
{
    int i, j, r;
    pthread_rwlock_init(&rwlock, NULL);   //初始化读写锁
    pthread_t* handles = (pthread_t*)malloc(NUM_THREADS * sizeof(pthread_t));// 创建对应的 Handle
    threadParam_t* param = (threadParam_t*)malloc(NUM_THREADS * sizeof(threadParam_t));// 创建对应的线程数据结构
    for (int t_id = 0; t_id < NUM_THREADS; t_id++)
    {
        param[t_id].t_id = t_id;
        param[t_id].matr = matr;
        param[t_id].n = n;
        param[t_id].m = m;
        param[t_id].l = l;
        pthread_create(&handles[t_id], 0, Static_threadFunc, &param[t_id]);
    }
    for (int t_id = 0; t_id < NUM_THREADS; t_id++)
        pthread_join(handles[t_id], 0);
}
void OpenMP(int n, int m, int l, int** matr)
{
    int flag;
#pragma omp parallel for shared(matr,n,m) private(flag)
    for (int id = 0;id < l;id++)
    {
        flag = 0;
        int* a = (int*)aligned_alloc(32, n * sizeof(int));
        int temp = shuru(a, n, m);
        for (int i = 0; i < n; i++)
        {
            if (flag == 1) continue;
            int bpos = i / 32, spos = i % 32;
            if (a[bpos] >> (31 - spos) != 0)
            {
                if (matr[n - 1 - i][bpos] >> (31 - spos) == 0)
                {
#pragma omp critical
                    for (int j = 0; j < m; j++)
                        matr[n - 1 - i][j] = a[j];
                    flag = 1;
                }
                else
                {
                    for (int j = 0; j < m; j++)
                        a[j] ^= matr[n - 1 - i][j];
                }
            }
        }
    }
}
void Parallel(int n, int m, int* a, int** matr)
{
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        int bpos = i / 32, spos = i % 32;
        if (a[bpos] >> (31 - spos) != 0)
        {
            if (matr[n - 1 - i][bpos] >> (31 - spos) == 0)
            {
                for (j = 0; j <= m - 4; j += 4)
                {
                    int32x4_t temp = vld1q_s32(a + j);
                    vst1q_s32(matr[n - 1 - i] + j, temp);
                }
                for (; j < m; j++)
                    matr[n - 1 - i][j] = a[j];
                return;
            }
            else
            {
                for (j = 0; j <= m - 4; j += 4)
                {
                    int32x4_t t1 = vld1q_s32(matr[n - 1 - i] + j);
                    int32x4_t t2 = vld1q_s32(a + j);
                    t2 = veorq_s32(t1, t2);
                    vst1q_s32(a + j, t2);
                }
                for (; j < m; j++)
                    a[j] ^= matr[n - 1 - i][j];
            }
        }
    }
}
void Serial(int n, int m, int* a, int** matr)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        int bpos = i / 32, spos = i % 32;
        if (a[bpos] >> (31 - spos) != 0)
        {
            if (matr[n - 1 - i][bpos] >> (31 - spos) == 0)
            {
                for (j = 0; j < m; j++)
                    matr[n - 1 - i][j] = a[j];
                return;
            }
            else
            {
                for (j = 0; j < m; j++)
                    a[j] ^= matr[n - 1 - i][j];
            }
        }
    }
}
int main()
{
    int n, m, u, l;
    srand(time(0));
    n = 200;
    while (n <= 10000)
    {
        u = rand() % (n / 2);
        while (u == 0) u = rand() % (n / 2);
        l = rand() % n + n / 2;
        if (n % 32 == 0)
            m = n / 32;
        else
            m = n / 32 + 1;
        int** ma = (int**)aligned_alloc(32, n * sizeof(int*));
        for (int i = 0; i < n; i++)
            ma[i] = (int*)aligned_alloc(32, m * sizeof(int));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                ma[i][j] = 0;
        for (int i = 0; i < u; i++)
        {
            int* a = (int*)aligned_alloc(32, n * sizeof(int));
            int temp = shuru(a, n, m);
            for (int j = 0; j < m; j++)
                ma[n - 1 - temp][j] = a[j];
        }

        struct timeval tstart, tend;
        double timeUsed1 = 0.0, timeUsed2 = 0.0, timeUsed3 = 0.0, timeUsed0 = 0.0, timeUsed4 = 0.0;

        for (int i = 0; i < l; i++)
        {
            gettimeofday(&tstart, NULL);
            int* a = (int*)aligned_alloc(32, n * sizeof(int));
            int temp = shuru(a, n, m);
            gettimeofday(&tend, NULL);
            timeUsed0 += 1000000 * (tend.tv_sec - tstart.tv_sec) + tend.tv_usec - tstart.tv_usec;

            int** mt = (int**)aligned_alloc(32, n * sizeof(int*));
            for (int i = 0; i < n; i++)
                mt[i] = (int*)aligned_alloc(32, m * sizeof(int));
            int* at = (int*)aligned_alloc(32, n * sizeof(int));

            for (int j = 0; j < n; j++)
                for (int k = 0; k < m; k++)
                    mt[j][k] = ma[j][k];
            for (int j = 0; j < n; j++)
                at[j] = a[j];
            gettimeofday(&tstart, NULL);
            Serial(n, m, at, mt);
            gettimeofday(&tend, NULL);
            timeUsed1 += 1000000 * (tend.tv_sec - tstart.tv_sec) + tend.tv_usec - tstart.tv_usec;

            for (int j = 0; j < n; j++)
                for (int k = 0; k < m; k++)
                    mt[j][k] = ma[j][k];
            for (int j = 0; j < n; j++)
                at[j] = a[j];
            gettimeofday(&tstart, NULL);
            Parallel(n, m, at, mt);
            gettimeofday(&tend, NULL);
            timeUsed2 += 1000000 * (tend.tv_sec - tstart.tv_sec) + tend.tv_usec - tstart.tv_usec;

            delete[]a;
            delete[]at;
            for (int i = 0;i < n;i++)
                delete[]mt[i];
            delete[]mt;
        }
        timeUsed0 /= 1000;
        timeUsed1 /= 1000;
        timeUsed2 /= 1000;

        gettimeofday(&tstart, NULL);
        Static_Pthread(n, m, l, ma);
        gettimeofday(&tend, NULL);
        timeUsed3 += 1000000 * (tend.tv_sec - tstart.tv_sec) + tend.tv_usec - tstart.tv_usec;
        timeUsed3 /= 1000;

        gettimeofday(&tstart, NULL);
        OpenMP(n, m, l, ma);
        gettimeofday(&tend, NULL);
        timeUsed4 += 1000000 * (tend.tv_sec - tstart.tv_sec) + tend.tv_usec - tstart.tv_usec;
        timeUsed4 /= 1000;

        cout << n << " " << timeUsed1 + timeUsed0 << " " << timeUsed2 + timeUsed0 << " " << timeUsed3 << " " << timeUsed4 << endl;
        if (n >= 3000) n += 1000;
        if (n >= 1000 && n < 3000) n += 500;
        if (n >= 500 && n < 1000) n += 100;
        if (n >= 200 && n < 500) n += 50;
    }
    pthread_rwlock_destroy(&rwlock);      //销毁读写锁
    return 0;
}