#include "LowerTriangularMatrix.h"
#include "Matrix.h"
#include <iostream>
#include <ctime>
#include <stdlib.h>

const size_t SIZE = 1440,
RANGE = 200;
const int START = -100;

inline void alalyzeNaive()
{
	Matrix<int> *naiveA = new Matrix<int>(SIZE),
		*naiveB = new Matrix<int>(SIZE);

	for (size_t i = 0; i < SIZE; i++)
		for (size_t j = 0; j <= i; j++)
		{
			int rnd = rand() % RANGE + START;
			naiveA->set(i, j, rnd);
			naiveA->set(j, i, rnd);
		}

	for (size_t i = 0; i < SIZE; i++)
		for (size_t j = 0; j <= i; j++)
			naiveB->set(i, j, rand() % RANGE + START);

	time_t start = clock();
	delete (*naiveA * naiveB);
	std::cout << "With naive matrix implementation: " << double(clock() - start) / CLOCKS_PER_SEC << std::endl << std::endl;
}

inline void analyze(size_t block_size)
{
	BlockMatrix<int> *blockA = new BlockMatrix<int>(SIZE / block_size, SIZE / block_size, block_size),
		*blockB = new BlockMatrix<int>(SIZE / block_size, SIZE / block_size, block_size);

	SymmetricMatrix<int> *symmA = new SymmetricMatrix<int>(SIZE, block_size);
	LowerTriangularMatrix<int> *lowerB = new LowerTriangularMatrix<int>(SIZE, block_size);

	// Заполнение случайными числами симметричной матрицы A
	for (size_t i = 0; i < SIZE; i++)
		for (size_t j = 0; j <= i; j++)
		{
			int rnd = rand() % RANGE + START;
			blockA->set(i, j, rnd);
			blockA->set(j, i, rnd);
			symmA->set(i, j, rnd);
		}

	// Заполнение случайными числами нижнетреугольной матрицы B
	for (size_t i = 0; i < SIZE; i++)
		for (size_t j = 0; j <= i; j++)
		{
			int rnd = rand() % RANGE + START;
			blockB->set(i, j, rnd);
			lowerB->set(i, j, rnd);
		}

	time_t start = clock();
	delete (*blockA * blockB);
	std::cout << "With block matrix implementation: " << double(clock() - start) / CLOCKS_PER_SEC << std::endl;

	start = clock();
	delete multiply(symmA, lowerB);
	std::cout << "With special matrix implementation: " << double(clock() - start) / CLOCKS_PER_SEC << std::endl;

	start = clock();
	delete multiplyWithOmp1(symmA, lowerB, 2);
	std::cout << "With special matrix implementation & 1st parallel multiplication: " << double(clock() - start) / CLOCKS_PER_SEC << std::endl;

	start = clock();
	delete multiplyWithOmp2(symmA, lowerB, 2);
	std::cout << "With special matrix implementation & 2nd parallel multiplication: " << double(clock() - start) / CLOCKS_PER_SEC << std::endl;
}

int main()
{
	std::cout << "Matrices dimension: " << SIZE << std::endl;
	srand(time(0));
	alalyzeNaive();
	unsigned int blockSizes[] = {6, 10, 15, 20, 24, 36, 40, 60, 72, 80, 96, 120, 144, 160, 180, 240, 360, 480, 720};

	//for (size_t bSize = 8; bSize < SIZE / 2; bSize *= 2)
	for (size_t i = 0; i < 19; i++)
	{
		std::cout << "Block size = " << blockSizes[i] << std::endl;
		analyze(blockSizes[i]);
		std::cout << std::endl;
	}	

	system("pause");
	return 0;
}