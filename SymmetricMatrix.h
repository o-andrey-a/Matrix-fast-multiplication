#pragma once

#include "Matrix.h"

template<typename T>
class SymmetricMatrix
{
private:
	T* matrix;
	size_t dim,
		blockSize,
		blockDim;
	

public:

	SymmetricMatrix()
	{
		dim = 0;
		blockSize = 0;
		blockDim = 0;
	}
	SymmetricMatrix(Matrix<T> m, size_t bSize)
	{
		if (m.getCols() != m.getRows())
			return;
		if (m.getCols() % bSize || m.getRows() % bSize)
			return;
		
		dim = m.getCols();
		blockSize = bSize;
		matrix = new T[dim * (dim + 1) / 2];
		
		for (size_t i = 0; i < dim; i++)
			for (size_t j = 0; j <= i; j++)
				set(i, j, m.get(i, j));
	}
	SymmetricMatrix(size_t dimension, size_t bSize)
	{
		if (m.cols % bSize || m.rows % bSize)
			return;

		dim = dimension;
		blockSize = bSize;
		matrix = new T[dim * (dim + 1) / 2];

		for (size_t i = 0; i < dim; i++)
			for (size_t j = 0; j <= i; j++)
				matrix.set(i, j, 0);
	}
	~SymmetricMatrix()
	{
	//	delete[] matrix;
	}

	T get(size_t row, size_t col)
	{
		if (col > dim || row > dim)
		{
			std::cout << "Index out of range, result is fictive" << std::endl;
			return 0;
		}

		if (col > row)
			return get(col, row);
		return matrix[col % blockSize + (row % blockSize) * blockSize +
			blockSize * blockSize *
			(col % blockSize + (row / blockSize) * (row / blockSize + 1) / 2)];
	}
	void set(size_t row, size_t col, T value)
	{
		if (col > dim || row > dim)
		{
			std::cout << "Index out of range" << std::endl;
			return;
		}

		if (col > row)
			get(col, row) = value;
		matrix[col % blockSize + (row % blockSize) * blockSize +
			blockSize * blockSize *
			(col % blockSize + (row / blockSize) * (row / blockSize + 1) / 2)] = value;
	}

	// Возвращает блок (row, col)
	Matrix<T> getBlock(size_t row, size_t col)
	{
		if (row > blockDim || col > blockDim)
			return Matrix<T>();

		Matrix<T> block(blockDim, blockDim);

		// Вычисление номера ячейки, с которой начинается нужный блок
		size_t cell = blockSize * blockSize *
			(col % blockSize + (row / blockSize) * (row / blockSize + 1) / 2);

		for (size_t i = 0; i < blockSize * blockSize; i++, cell++)
			block.matrix[i] = matrix[cell];

		return block;
	}
	void setBlock(size_t row, size_t col, Matrix<T> m)
	{
		if (m.getCols() > cols - col ||
			m.getRows() > rows - row)
			return;

		size_t rowsStart = row * blockSize,
			rowsEnd = rowsStart + m.getRows(),
			colsStart = col * blockSize,
			colsEnd = colsStart + m.getCols();

		for (size_t i = rowsStart; i < rowsEnd; i++)
			for (size_t j = colsStart; j < colsEnd; j++)
				this->set(i, j, m.get(i - rowsStart, j - colsStart));
	}

	BlockMatrix<T> multiply(BlockMatrix<T> that)
	{
		if (this->blockCols != that.blockRows ||
			this->blockSize != that.blockSize)
			return BlockMatrix<T>();

		BlockMatrix<T> product(this->blockRows, that.blockCols, this->blockSize);

		for (size_t i = 0; i < this->blockRows; i++)
			for (size_t j = 0; j < that.blockCols; j++)
			{
				Matrix<T> accumulator(blockSize, blockSize);
				for (size_t k = 0; k < this->blockCols; k++)
					accumulator += this->getBlock(i, k) * that.getBlock(k, j);
				product.setBlock(i, j, accumulator);
			}

		return product;
	}
	BlockMatrix<T> multiply(SymmetricMatrix<T> that)
	{
		if (this->dim != that.dim)
			return BlockMatrix<T>();

		BlockMatrix<T> product(dim, dim, blockSize);

		for (size_t i = 0; i < this->dim; i++)
			for (size_t j = 0; j < that.dim; j++)
			{
				Matrix<T> accumulator(blockSize, blockSize);
				for (size_t k = 0; k < dim; k++)
					accumulator += this->getBlock(i, k) * that.getBlock(k, j);
				product.setBlock(i, j, accumulator);
			}

		return product;
	}

	void operator = (Matrix<T> m)
	{
		if (this->dim != m.cols ||
			this->dim != m.rows)
			//	Что-то сделать
			return;

		for (size_t i = 0; i < dim; i++)
			for (size_t j = 0; j <= i; j++)
				this->set(i, j, m.get(i, j));
	}
	void operator = (BlockMatrix<T> m)
	{
		if (this->dim != m.cols ||
			this->dim != m.rows)
			//	Что-то сделать
			return;

		for (size_t i = 0; i < dim; i++)
			for (size_t j = 0; j <= i; j++)
				this->set(i, j, m.get(i, j));
	}

	void print()
	{
		for (size_t i = 0; i < rows; i++)
		{
			for (size_t j = 0; j < cols; j++)
				std::cout << get(i, j) << "	";
			std::cout << "\n";
		}
		std::cout << std::endl;
	}
};

