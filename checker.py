import os
import sys
import random
import numpy as np
import time

n = int(sys.argv[1])

test_correctness = False

def gen_input(n):
    with open('input_rand.txt', 'w') as f:
        for i in range(n):
            for j in range(n):
                print(4*random.random()-2, end=' ', file=f) # in range [-2,2)
            print('', file=f)

def test(cond, msg):
    if cond:
        print(f'SUCCESS: {msg}')
    else:
        print(f'FAIL: {msg}')

if __name__=='__main__':
    gen_input(n)
    print(f'Random matrix of {n}X{n} generated.')
    os.system('./compile.sh')
    print(f'Running strategy 0....')
    start = time.time()
    os.system(f'bash run.sh {n} input_rand.txt 1 0')
    seq = time.time() - start
    for n_threads in 2,4,8,16:
        print(f'--- n_threads={n_threads}')
        for strategy in [1,2,3,4]:
            print(f'Running strategy {strategy}....')
            start = time.time()
            os.system(f'bash run.sh {n} input_rand.txt {n_threads} {strategy}')
            end = time.time()
            print(f'Time: {end-start}s, Speedup: {seq/(end-start)}')
            if test_correctness:
                A = np.loadtxt('input_rand.txt')
                L = np.loadtxt(f'output_L_{strategy}_{n_threads}.txt')
                U = np.loadtxt(f'output_U_{strategy}_{n_threads}.txt')
                test(np.allclose(L, np.tril(L)), 'L lower triangular')
                test(np.allclose(U, np.triu(U)), 'U upper triangular')
                test(abs(np.linalg.det(U) - 1) < 1e-3, 'Determinant of U less than 1')
                a_dash=np.matmul(L,U)
                test(np.allclose(a_dash, A, rtol=0, atol=1e-3), 'L*U identical to A')

        
