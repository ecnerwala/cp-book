#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <cstdio>
#include <cstring>
#include <string_view>
#include <type_traits>

namespace wala {

// Buffered writer with fast integer formatting, for problems whose output would otherwise dominate the runtime.
// flush() hands the buffer to the FILE* (not to the OS); it flushes on destruction, and anything else written
// through the same FILE* must be ordered with explicit flush() calls.
class FastWriter {
	static constexpr int BUF_SIZE = 1 << 16;
	static constexpr int MAX_ITEM = 32;
	FILE* file;
	char buf[BUF_SIZE + MAX_ITEM];
	int pos = 0;

	void reserve() {
		if (pos > BUF_SIZE) flush();
	}

public:
	explicit FastWriter(FILE* file_ = stdout) : file(file_) {}
	FastWriter(const FastWriter&) = delete;
	FastWriter& operator=(const FastWriter&) = delete;
	~FastWriter() { flush(); }

	void flush() {
		fwrite(buf, 1, pos, file);
		pos = 0;
	}

	FastWriter& operator<<(char c) {
		reserve();
		buf[pos++] = c;
		return *this;
	}

	FastWriter& operator<<(std::string_view s) {
		for (size_t i = 0; i < s.size(); i += MAX_ITEM) {
			reserve();
			int len = int(std::min(s.size() - i, size_t(MAX_ITEM)));
			memcpy(buf + pos, s.data() + i, len);
			pos += len;
		}
		return *this;
	}

	template <std::integral T> requires (!std::same_as<T, char> && !std::same_as<T, bool>)
	FastWriter& operator<<(T v) {
		static constexpr auto DIGITS = [] {
			std::array<char, 200> d{};
			for (int i = 0; i < 100; i++) d[2*i] = char('0' + i/10), d[2*i+1] = char('0' + i%10);
			return d;
		}();
		reserve();
		using U = std::make_unsigned_t<T>;
		U u = U(v);
		if constexpr (std::is_signed_v<T>) {
			if (v < 0) buf[pos++] = '-', u = U(0) - u;
		}
		char tmp[MAX_ITEM];
		int len = 0;
		for (; u >= 100; u /= 100) {
			len += 2;
			memcpy(tmp + MAX_ITEM - len, &DIGITS[2 * int(u % 100)], 2);
		}
		if (u >= 10) {
			len += 2;
			memcpy(tmp + MAX_ITEM - len, &DIGITS[2 * int(u)], 2);
		} else {
			tmp[MAX_ITEM - ++len] = char('0' + u);
		}
		memcpy(buf + pos, tmp + MAX_ITEM - len, len);
		pos += len;
		return *this;
	}
};

} // namespace wala
