Main function: decoptL1L1.m

This file is a modification of the original decoptSolver to make the algorithm
work for L1-L1 minimization:

	minimize    ||x||_1 + beta*||x - w||_1
        x
	subject to  A*x = b

where A, b, beta, and w are given.

The rest of the code was provided by Quoc Tran-Dinh and Volkan Cevher both at
Laboratory for Information and Inference Systems (LIONS) ad EPFL, Switzerland.

