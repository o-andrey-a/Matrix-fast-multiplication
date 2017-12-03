#pragma once
#include <ctime>
#include <iostream>
#include <omp.h>
#include <stdlib.h>

template<typename T>
class BlockMatrix;

template<typename T>
class SymmetricMatrix;

template<typename T>
class Matrix
{
private:
	T* matrix;
	size_t len,
		rows,
		cols;

	T operator [] (size_t k)
	{
		return matrix[k];
	}

public:	

	size_t getRows()
	{
		return rows;
	}

	size_t getCols()
	{
		return cols;
	}
	
	Matrix()
	{
		rows = 0;
		cols = 0;
		len = 0;
	}

	Matrix(size_t a)
	{
		rows = a;
		cols = a;
		len = rows * cols;
		matrix = new T[len];
		for (size_t i = 0; i < len; i++)
			matrix[i] = 0;
	}

	Matrix(size_t a, size_t b)
	{
		rows = a;
		cols = b;
		len = rows * cols;
		matrix = new T[len];
		for (size_t i = 0; i < len; i++)
			matrix[i] = 0;
	}

	~Matrix()
	{
		delete[] matrix;
	}

	void set(size_t row, size_t col, T value)
	{
		matrix[row * cols + col] = value;
	}

	T get(size_t row, size_t col)
	{
		return matrix[row * cols + col];
	}

	// Возвращает матрицу размера size x size,
	// для которой первый элемент расположен в ячейке (row, col)
	Matrix<T>* getBlock(size_t row, size_t col, size_t size)
	{
		Matrix<T>* block = new Matrix<T>(size);
		for (size_t i = 0; i < size; i++)
			for (size_t j = 0; j < size; j++)
				block->set(i, j, this->get(row + i, col + j));
		
		return block;
	}

	void operator = (Matrix<T>* m)
	{
		if (m->len != this->len)
		{
			len = m->len;
			rows = m->rows;
			cols = m->cols;
			matrix = new T[len];
		}

		for (size_t i = 0; i < len; i++)
			this->matrix[i] = m->matrix[i];
	}

	Matrix<T>* operator + (Matrix<T>* that)
	{
		if (this->cols != that->getCols() ||
			this->rows != that->getRows())
			return (new Matrix<T>());

		Matrix<T>* sum = new Matrix<T>(rows, cols);
		for (size_t i = 0; i < len; i++)
			sum->matrix[i] = this->matrix[i] + that->matrix[i];

		return sum;
	}

	void operator += (Matrix<T>* that)
	{
		if (this->cols != that->getCols() ||
			this->rows != that->getCols())
		{
			std::cout << "Matrix dimensions must agree" << std::endl;
			return;
		}

		for (size_t i = 0; i < len; i++)
			this->matrix[i] = this->matrix[i] + that->matrix[i];
	}

	Matrix<T>* operator * (Matrix<T>* that)
	{
		if (this->cols != that->getRows())
			return (new Matrix<T>());

		Matrix<T>* product = new Matrix(rows, that->cols);
		for (size_t i = 0; i < this->rows; i++)
			for (size_t j = 0; j < that->cols; j++)
			{
				T val = 0;
				for (size_t k = 0; k < this->cols; k++)
					val += this->get(i, k) * that->get(k, j);
				product->set(i, j, val);
			}
		return product;
		
	}

	Matrix<T>* multiplyWithOmp(Matrix<T>* that, size_t numThreads)
	{
		if (this->cols != that->getRows())
			return (new Matrix<T>());

		Matrix<T>* product = new Matrix(rows, that->cols);

#pragma omp parallel for num_threads(numThreads)
		for (int i = 0; i < this->rows; i++)
			for (size_t j = 0; j < that->cols; j++)
			{
				T val = 0;
				for (size_t k = 0; k < this->cols; k++)
					val += this->get(i, k) * that->get(k, j);
				product->set(i, j, val);
			}
		return product;
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

	// Генерирует матрицу случайных чисел в диапазоне от start до start + range
	void generateRandom(T start, T range)
	{
		for (size_t i = 0; i < len; i++)
			matrix[i] = rand() % range + start;
	}

	friend BlockMatrix<T>;
	friend SymmetricMatrix<T>;

	bool operator == (Matrix<T>* m)
	{
		for (size_t i = 0; i < rows; i++)
			for (size_t j = 0; j < cols; j++)
				if (this->get(i, j) != m->get(i, j))
					return false;

		return true;
	}

};







template<typename T>
class BlockMatrix
{

private:
	T* matrix;

public:
	size_t blockSize,
		blockRows,
		blockCols,
		rows,
		cols;

	BlockMatrix()
	{
		blockSize = 0;
		blockRows = 0;
		blockCols = 0;
	}

	BlockMatrix(Matrix<T>* m, size_t blockVolume)
	{
		if (m->cols % blockVolume || m->rows % blockVolume)
		{
			std::cout << "Can't split matrix into blocks of this dimension" << std::endl;
			return;
		}

		cols = m->cols;
		rows = m->rows;
		blockSize = blockVolume;
		blockCols = cols / blockVolume;
		blockRows = rows / blockVolume;
		matrix = new T[rows * cols];

		for (size_t i = 0; i < blockRows; i++)
			for (size_t j = 0; j < blockCols; j++)
			{
				Matrix<T>* deleteMe = m->getBlock(i * blockVolume, j * blockVolume, blockVolume);
				this->setBlock(i, j, deleteMe);
				delete deleteMe;
			}
	}

	BlockMatrix(size_t bRows, size_t bCols, size_t bSize)
	{
		blockRows = bRows;
		blockCols = bCols;
		blockSize = bSize;
		rows = bRows * bSize;
		cols = bCols * bSize;
		matrix = new T[rows * cols];

		for (size_t i = 0; i < rows * cols; i++)
			matrix[i] = 0;
	}

	~BlockMatrix()
	{
		delete[] matrix;
	}

	T get(size_t row, size_t col)
	{
		// Матрица хранится по блочным строкам
		return matrix[blockSize * (row / blockSize * blockSize * blockCols
			+ col / blockSize * blockSize + row % blockSize)
			+ col % blockSize];
	}

	void set(size_t row, size_t col, T value)
	{
		matrix[blockSize * (row / blockSize * blockSize * blockCols
			+ col / blockSize * blockSize + row % blockSize)
			+ col % blockSize] = value;
	}

	// Возвращает блок (row, col)
	Matrix<T>* getBlock(size_t row, size_t col)
	{
		if (row > blockRows || col > blockCols)
			return (new Matrix<T>());

		Matrix<T>* block = new Matrix<T>(blockSize, blockSize);

		// Вычисление номера ячейки, с которой начинается нужный блок
		size_t cell = blockSize * blockSize * (row * blockCols + col);

		for (size_t i = 0; i < blockSize * blockSize; i++, cell++)
			block->matrix[i] = matrix[cell];

		return block;
	}

	void setBlock(size_t row, size_t col, Matrix<T>* m)
	{
		if (m->getCols() > cols - col ||
			m->getRows() > rows - row)
			return;

		size_t rowsStart = row * blockSize,
			rowsEnd = rowsStart + m->getRows(),
			colsStart = col * blockSize,
			colsEnd = colsStart + m->getCols();

		for (size_t i = rowsStart; i < rowsEnd; i++)
			for (size_t j = colsStart; j < colsEnd; j++)
				this->set(i, j, m->get(i - rowsStart, j - colsStart));
	}

	BlockMatrix<T>* operator * (BlockMatrix<T>* that)
	{
		if (this->blockCols != that->blockRows ||
			this->blockSize != that->blockSize)
			return (new BlockMatrix<T>());

		BlockMatrix<T>* product = new BlockMatrix<T>(this->blockRows, that->blockCols, this->blockSize);

		for (size_t i = 0; i < this->blockRows; i++)
			for (size_t j = 0; j < that->blockCols; j++)
			{
				Matrix<T> accumulator(blockSize, blockSize);
				Matrix<T>* deleteMe = &accumulator;
				for (size_t k = 0; k < this->blockCols; k++)
				{
					Matrix<T> *factor1 = this->getBlock(i, k),
						*factor2 = that->getBlock(k, j),
						*res = (*factor1 * factor2);
					
					accumulator += res;

					delete factor1;
					delete factor2;
					delete res;
				}
				product->setBlock(i, j, deleteMe);
				
			}

		return product;
	}

	void operator = (Matrix<T>* m)
	{
		if (this->cols != m->cols ||
			this->rows != m->rows)
			//	Что-то сделать
			return;

		for (size_t i = 0; i < rows; i++)
			for (size_t j = 0; j < cols; j++)
				this->set(i, j, m->get(i, j));
	}

	void operator = (BlockMatrix<T>* m)
	{
		if (this->cols != m->cols ||
			this->rows != m->rows)
			//	Что-то сделать
			return;

		for (size_t i = 0; i < rows * cols; i++)
			matrix[i] = m->matrix[i];
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

	void generateRandom(T start, T range)
	{
		for (size_t i = 0; i < rows * cols; i++)
			matrix[i] = rand() % range + start;
	}

	bool operator == (Matrix<T>* m)
	{
		for (size_t i = 0; i < rows; i++)
			for (size_t j = 0; j < cols; j++)
				if (this->get(i, j) != m->get(i, j))
					return false;
		return true;
	}

	bool operator == (BlockMatrix<T>* m)
	{
		for (size_t i = 0; i < rows; i++)
			for (size_t j = 0; j < cols; j++)
				if (this->get(i, j) != m->get(i, j))
					return false;
		return true;
	}

};







template<typename T>
class SymmetricMatrix
{
private:
	T* matrix;

public:
	size_t dim,
		blockSize,
		blockDim,
		len;

	SymmetricMatrix()
	{
		dim = 0;
		blockSize = 0;
		blockDim = 0;
		len = 0;
	}

	SymmetricMatrix(Matrix<T> m, size_t bSize)
	{
		if (m->getCols() != m->getRows())
			return;
		if (m->getCols() % bSize || m->getRows() % bSize)
			return;

		dim = m->getCols();
		blockSize = bSize;
		blockDim = dim / blockSize;
		len = dim * (dim + 1) / 2;
		matrix = new T[len];

		for (size_t i = 0; i < dim; i++)
			for (size_t j = 0; j <= i; j++)
				set(i, j, m->get(i, j));
	}

	SymmetricMatrix(size_t dimension, size_t bSize)
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

	~SymmetricMatrix()
	{
		delete[] matrix;
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

		size_t elementsAbove,
			elementsToTheLeft = blockSize * blockSize * (col / blockSize),
			posInBlock,
			rowInBlock = row % blockSize;

		posInBlock = col % blockSize + (col < row - rowInBlock ?
			rowInBlock * blockSize : rowInBlock * (rowInBlock + 1) / 2);

		elementsAbove = (row - rowInBlock) * (row - rowInBlock + 1) / 2;

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
			set(col, row, value);
			return;
		}

		size_t elementsAbove,
			elementsToTheLeft = blockSize * blockSize * (col / blockSize),
			posInBlock,
			rowInBlock = row % blockSize;

		posInBlock = col % blockSize + (col < row - rowInBlock ?
			rowInBlock * blockSize : rowInBlock * (rowInBlock + 1) / 2);

		elementsAbove = (row - rowInBlock) * (row - rowInBlock + 1) / 2;

		matrix[elementsAbove + elementsToTheLeft + posInBlock] = value;
	}

	// Возвращает блок (row, col)
	Matrix<T>* getBlock(size_t row, size_t col)
	{
		if (row > blockDim || col > blockDim)
			return (new Matrix<T>());

		Matrix<T> *block = new Matrix<T>(blockSize, blockSize);

		if (col > row)
		{
			for (size_t i = 0; i < blockSize; i++)
				for (size_t j = 0; j < blockSize; j++)
					block->set(i, j, this->get(col * blockSize + j, row * blockSize + i));
		}
		else
		{
			for (size_t i = 0; i < blockSize; i++)
				for (size_t j = 0; j < blockSize; j++)
					block->set(i, j, this->get(row * blockSize + i, col * blockSize + j));
		}
		return block;
	}

	void setBlock(size_t row, size_t col, Matrix<T> m)
	{
		if (m->getCols() > dim - col ||
			m->getRows() > dim - row)
			return;

		size_t rowsStart = row * blockSize,
			rowsEnd = rowsStart + m->getRows(),
			colsStart = col * blockSize,
			colsEnd = colsStart + m->getCols();

		for (size_t i = rowsStart; i < rowsEnd; i++)
			for (size_t j = colsStart; j < colsEnd; j++)
				this->set(i, j, m->get(i - rowsStart, j - colsStart));
	}

	BlockMatrix<T>* multiply(BlockMatrix<T>* that)
	{
		if (this->blockDim != that->blockRows ||
			this->blockSize != that->blockSize)
			return (new BlockMatrix<T>());

		BlockMatrix<T>* product = new BlockMatrix<T>(this->blockDim, that->blockCols, this->blockSize);

		for (size_t i = 0; i < this->blockDim; i++)
			for (size_t j = 0; j < that->blockCols; j++)
			{
				Matrix<T> accumulator(blockSize, blockSize);
				Matrix<T>* deleteMe = &accumulator;
				for (size_t k = 0; k < this->blockDim; k++)
				{
					Matrix<T> *factor1 = this->getBlock(i, k),
						*factor2 = that->getBlock(k, j),
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

	BlockMatrix<T>* multiply(SymmetricMatrix<T>* that)
	{
		if (this->dim != that->dim)
			return (new BlockMatrix<T>());

		BlockMatrix<T> product = new BlockMatrix<T>(blockDim, blockDim, blockSize);

		for (size_t i = 0; i < blockDim; i++)
			for (size_t j = 0; j < blockDim; j++)
			{
				Matrix<T> accumulator(blockSize, blockSize);
				Matrix<T>* deleteMe = &accumulator;
				for (size_t k = 0; k < this->blockDim; k++)
				{
					Matrix<T> *factor1 = this->getBlock(i, k),
						*factor2 = that->getBlock(k, j),
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

	void operator = (Matrix<T>* m)
	{
		if (this->dim != m->cols ||
			this->dim != m->rows)
			//	Что-то сделать
			return;

		for (size_t i = 0; i < dim; i++)
			for (size_t j = 0; j <= i; j++)
				this->set(i, j, m->get(i, j));
	}

	void operator = (BlockMatrix<T>* m)
	{
		if (this->dim != m->cols ||
			this->dim != m->rows)
			//	Что-то сделать
			return;

		for (size_t i = 0; i < dim; i++)
			for (size_t j = 0; j <= i; j++)
				this->set(i, j, m->get(i, j));
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

	void generateRandom(T start, T range)
	{
		for (size_t i = 0; i < len; i++)
			matrix[i] = rand() % range + start;
	}

	bool operator == (Matrix<T>* m)
	{
		for (size_t i = 0; i < dim; i++)
			for (size_t j = 0; j < dim; j++)
				if (this->get(i, j) != m->get(i, j))
					return false;
		return true;
	}


};