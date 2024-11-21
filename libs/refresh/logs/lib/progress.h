#pragma once

#include <cinttypes>
#include <string>
#include <cmath>
#include <charconv>

namespace refresh
{
	class progress_state
	{
		enum class type_t {none, number, of_total, percent};

		type_t type{ type_t::none };

		uint64_t counter{ 0 };
		uint64_t total{ 0 };

		std::string prefix;
		std::string separator;
		std::string suffix;

		int precision{ -1 };
		double mult{ 1 };

		std::string message;
		bool message_checked{ false };

		void adjust_precision()
		{
			if (total == 0)
			{
				total = 1;				// Should be error message here?
			}

			if (precision < 0)
			{
				if (total <= 100)
					precision = 0;
				else if (total <= 10000)
					precision = 1;
				else if (total <= 1000000)
					precision = 2;
				else
					precision = 3;
			}
			else if (precision > 6)
				precision = 6;
				
			mult = pow(10, precision);
		}

		void build_message()
		{
			std::string s;
			char buffer[16];

			if (type == type_t::number)
				s = std::to_string(counter);
			else if (type == type_t::of_total)
				s = prefix + std::to_string(counter) + separator + std::to_string(total) + suffix;
			else if (type == type_t::percent)
			{
				auto r = std::to_chars(buffer, buffer + 16, 100.0 * counter / total, std::chars_format::fixed, precision);
				
				if (r.ec != std::errc())
					s = "";
				else
				{
					s = prefix;
					s.append(buffer, r.ptr - buffer);
					s += suffix;
				}
			}

			if (s != message)
			{
				message = move(s);
				message_checked = false;
			}
		}

	public:
		progress_state() :
			type(type_t::number)
		{}

		progress_state(uint64_t total, const std::string &prefix, const std::string &separator, const std::string &suffix) :
			type(type_t::of_total),
			total(total),
			prefix(prefix),
			separator(separator),
			suffix(suffix)
		{}

		progress_state(uint64_t total, const std::string& prefix, const std::string& suffix, int precision) :
			type(type_t::percent),
			total(total),
			prefix(prefix),
			suffix(suffix),
			precision(precision)
		{
			adjust_precision();
		}

		const std::string& str()
		{
			message_checked = true;
			return message;
		}

		bool increment(uint64_t n)
		{
			counter += n;
			build_message();

			return !message_checked;
		}

		bool was_checked() const
		{
			return message_checked;
		}
	};
}