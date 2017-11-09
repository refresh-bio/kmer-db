#pragma once
#include <vector>

template <class T>
class Array
{

public:
	size_t size() const { return width; }
	const std::vector<T>& getData() const { return v; }
	std::vector<T>& getData() { return v; }

	Array(int size) : width(size), height(size), v(size * size) {}
	Array() : width(0), height(0) {}
	Array(const Array<T>& ref) : width(ref.width), height(ref.height), v(ref.v) {}
	Array(Array&& rhs) : width(rhs.width), height(rhs.height), v(std::move(rhs.v)) {}

	Array<T>& operator=(const Array<T>& ref) = default;
	Array<T> operator-(const Array<T>& b) {
		Array<T> out;
		if (this->size() == b.size()) {
			out.resize(this->size(), this->size(), 0);

			for (int i = 0; i < this->size() * this->size(); ++i) {
				out.getData()[i] = this->getData()[i] - b.getData()[i];
			}
		}
		return out;
	}

	void resize(int width, int height) {
		this->width = width;
		this->height = height;
		v.resize(width * height);
	}

	void resize(int width, int height, const T& value) {
		this->width = width;
		this->height = height;
		v.resize(width * height, value);
	}

	void clear() {
		v.clear();
	}

	T* operator[](const int row) { return v.data() + row * width; }
	const T* operator[](const int row) const { return v.data() + row * width; }

	void save(std::ofstream & file) {
		for (int i = 0; i < size(); ++i) {
			for (int j = 0; j < size(); ++j) {
				file  << (*this)[i][j] << ',';
			}
			file << endl;
		}
	}


protected:
	int width;
	int height;
	std::vector<T> v;
};