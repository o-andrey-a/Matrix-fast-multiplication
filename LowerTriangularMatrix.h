#pragma once
#include "Matrix.h"
#include <omp.h>

template<typename T>
class LowerTriangularMatrix
{
private:
	T* matrix;

public:
	size_t dim,
		blockSize,
		blockDim,
		len;

	LowerTriangularMatrix()
	{
		dim = 0;
		blockSize = 0;
		blockDim = 0;
		len = 0;
	}

	LowerTriangularMatrix(size_t dimension, size_t bSize)
	{
		if (dimension % bSize)
			return;

		dim = dimension;
		blockSize = bSize;
		blockDim = dim / blockSize;
		len = dim * (dim + 1) / 2;
		matrix = new T[len];

		for (size_t i = 0; i < len; i++)
			matrix[i] = 0;
	}

	T get(size_t row, size_t col)
	{
		if (col > dim || row > dim)
		{
			std::cout << "Index out of range, result is fictive" << std::endl;
			return 0;
		}

		if (col > row)
			return 0;

		size_t elementsAbove,
			elementsToTheLeft,
			posInBlock,
			rowInBlock = row % blockSize,
			colInBlock = col % blockSize,
			colsToTheLeft = blockSize * (col / blockSize);

		elementsToTheLeft = (2 * dim - colsToTheLeft + 1) * colsToTheLeft / 2;

		if (col < row - rowInBlock)
		{
			posInBlock = colInBlock + rowInBlock * blockSize;
			elementsAbove = blockSize * (blockSize + 1) / 2 +

				(row - col - rowInBlock - (blockSize - colInBlock)) / blockSize
				* blockSize * blockSize;
		}
		else
		{
			posInBlock = colInBlock + rowInBlock * (rowInBlock + 1) / 2;
			elementsAbove = 0;
		}

		return matrix[elementsAbove + elementsToTheLeft + posInBlock];
	}

	void set(size_t row, size_t col, T value)
	{
		if (col > dim || row > dim)
		{
			std::cout << "Index out of range" << std::endl;
			return;
		}

		if (col > row)
		{
			std::cout << "You're trying to modify the upper triangle of the lower triangular matrix" << std::endl;
			return;
		}

		size_t elementsAbove,
			elementsToTheLeft,
			posInBlock,
			rowInBlock = row % blockSize,
			colInBlock = col % blockSize,
			colsToTheLeft = blockSize * (col / blockSize);

		elementsToTheLeft = (2 * dim - colsToTheLeft + 1) * colsToTheLeft / 2;

		if (col < row - rowInBlock)
		{
			posInBlock = colInBlock + rowInBlock * blockSize;
			elementsAbove = blockSize * (blockSize + 1) / 2 +

				(row - col - rowInBlock - (blockSize - colInBlock)) / blockSize
				* blockSize * blockSize;
		}
		else
		{
			posInBlock = colInBlock + rowInBlock * (rowInBlock + 1) / 2;
			elementsAbove = 0;
		}

		matrix[elementsAbove + elementsToTheLeft + posInBlock] = value;
	}

	Matrix<T>* getBlock(size_t row, size_t col)
	{
		if (row >= blockDim || col >= blockDim)
			return (new Matrix<T>());

		Matrix<T> *block = new Matrix<T>(blockSize, blockSize);

		if (col <= row)
		{
			for (size_t i = 0; i < blockSize; i++)
				for (size_t j = 0; j < blockSize; j++)
					block->set(i, j, this->get(row * blockSize + i, col * blockSize + j));
		}
		return block;
	}

	BlockMatrix<T>* operator * (BlockMatrix<T>* that)
	{
		if (this->blockDim != that->blockRows)
			return (new BlockMatrix<T>());

		BlockMatrix<T>* product = new BlockMatrix<T>(this->blockDim, that->blockCols, blockSize);

		for (size_t i = 0; i < blockDim; i++)
			for (size_t j = 0; j < that->blockCols; j++)
			{
				Matrix<T>* accumulator = new Matrix<T>(blockSize, blockSize);
				for (size_t k = 0; k <= i; k++)
					accumulator += this->getBlock(i, k) * that->getBlock(k, j);
				product->setBlock(i, j, accumulator);
			}

		return product;
	}

	bool operator == (Matrix<T>* m)
	{
		for (size_t i = 0; i < dim; i++)
			for (size_t j = 0; j < dim; j++)
				if (this->get(i, j) != m->get(i, j))
					return false;
		return true;
	}

	void print()
	{
		for (size_t i = 0; i < dim; i++)
		{
			for (size_t j = 0; j < dim; j++)
				std::cout << get(i, j) << "	";
			std::cout << "\n";
		}
		std::cout << std::endl;
	}

};

template<typename T>
BlockMatrix<T>* multiply(SymmetricMatrix<T>* m1, LowerTriangularMatrix<T>* m2)
{
	if (m1->blockDim != m2->blockDim)
		return (new BlockMatrix<T>());

	BlockMatrix<T>* product = new BlockMatrix<T>(m1->blockDim, m2->blockDim, m1->blockSize);

	for (size_t i = 0; i < m1->blockDim; i++)
		for (size_t j = 0; j < m2->blockDim; j++)
		{
			Matrix<T> accumulator(m1->blockSize, m1->blockSize);
			Matrix<T>* deleteMe = &accumulator;
			for (size_t k = 0; k <= i; k++)
			{
				Matrix<T> *factor1 = m1->getBlock(i, k),
					*factor2 = m2->getBlock(k, j),
					*res = (*factor1 * factor2);

				product->setBlock(i, j, res);

				delete factor1;
				delete factor2;
				delete res;
			}
			product->setBlock(i, j, deleteMe);
		}
	return product;
}

template<typename T>
BlockMatrix<T>* multiplyWithOmp1(SymmetricMatrix<T>* m1, LowerTriangularMatrix<T>* m2, size_t numThreads)
{
	if (m1->blockDim != m2->blockDim)
		return (new BlockMatrix<T>());

	BlockMatrix<T>* product = new BlockMatrix<T>(m1->blockDim, m2->blockDim, m1->blockSize);

#pragma omp parallel for num_threads(numThreads)
	for (int i = 0; i < m1->blockDim; i++)
		for (size_t j = 0; j < m2->blockDim; j++)
		{
			Matrix<T> accumulator(m1->blockSize, m1->blockSize);
			Matrix<T>* deleteMe = &accumulator;
			for (size_t k = 0; k <= i; k++)
			{
				Matrix<T> *factor1 = m1->getBlock(i, k),
					*factor2 = m2->getBlock(k, j),
					*res = (*factor1 * factor2);

				product->setBlock(i, j, res);

				delete factor1;
				delete factor2;
				delete res;
			}
			product->setBlock(i, j, deleteMe);
		}

	return product;
}

template<typename T>
BlockMatrix<T>* multiplyWithOmp2(SymmetricMatrix<T>* m1, LowerTriangularMatrix<T>* m2, size_t numThreads)
{
	if (m1->blockDim != m2->blockDim)
		return (new BlockMatrix<T>());

	BlockMatrix<T>* product = new BlockMatrix<T>(m1->blockDim, m2->blockDim, m1->blockSize);

	for (int i = 0; i < m1->blockDim; i++)
		for (size_t j = 0; j < m2->blockDim; j++)
		{
			Matrix<T> accumulator(m1->blockSize, m1->blockSize);
			Matrix<T>* deleteMe = &accumulator;
			for (size_t k = 0; k <= i; k++)
			{
				Matrix<T> *factor1 = m1->getBlock(i, k),
					*factor2 = m2->getBlock(k, j),
					*res = factor1->multiplyWithOmp(factor2, numThreads);

				product->setBlock(i, j, res);

				delete factor1;
				delete factor2;
				delete res;
			}
			product->setBlock(i, j, deleteMe);
		}

	return product;
}