#pragma once
#include <stdexcept>
#include <string>
namespace raduls
{
	namespace exceptions
	{
		class InputNotAlignedException : public std::logic_error
		{
		public:
			InputNotAlignedException(uint64_t alignment) :std::logic_error("Input array is not properly aligned to " + std::to_string(alignment))
			{}
		};

		class TempNotAlignedException : public std::logic_error
		{
		public:
			TempNotAlignedException(uint64_t alignment) :std::logic_error("Temporary array is not properly aligned to " + std::to_string(alignment))
			{}
		};

		class RecSizeNotMultipleOf8Exception : public std::logic_error
		{
		public:
			RecSizeNotMultipleOf8Exception() : std::logic_error("Rec size must be a multiple of 8")
			{}
		};

		class KeySizeGreaterThanRecSizeException : public std::logic_error
		{
		public:
			KeySizeGreaterThanRecSizeException() : std::logic_error("Key size cannot be greater than rec size")
			{}
		};

		class UsupportedRecSizeException : public std::logic_error
		{
		public:
			UsupportedRecSizeException() : std::logic_error("This rec size is not supported. Try to extend MAX_REC_SIZE_IN_BYTES")
			{}
		};
	}
}
