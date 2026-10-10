#include <catch2/catch_test_macros.hpp>

#include "tensor.hpp"

#include <algorithm>
#include <functional>
#include <iterator>
#include <memory>
#include <ranges>
#include <span>
#include <string>
#include <utility>
#include <vector>

using wala::vec;
using wala::bounded_vec;
using wala::bounded_stack;
using wala::with_capacity;
using wala::uninit;

template <typename C> constexpr bool container_ok =
	std::contiguous_iterator<typename C::iterator>
	&& std::contiguous_iterator<typename C::const_iterator>
	&& std::convertible_to<typename C::iterator, typename C::const_iterator>
	&& !std::convertible_to<typename C::const_iterator, typename C::iterator>
	&& std::ranges::contiguous_range<C>
	&& std::ranges::contiguous_range<const C>
	&& std::ranges::sized_range<C>
	&& std::same_as<std::ranges::range_size_t<C>, int>
	&& std::same_as<std::ranges::range_difference_t<C>, int>
	&& std::same_as<std::ranges::range_value_t<const C>, typename C::value_type>
	&& std::same_as<std::ranges::range_reference_t<const C>, const typename C::value_type&>
	&& std::is_nothrow_move_constructible_v<C>
	&& !std::is_copy_constructible_v<C>;
static_assert(container_ok<vec<int>>);
static_assert(container_ok<vec<std::string>>);
static_assert(container_ok<bounded_vec<int>>);
static_assert(container_ok<bounded_vec<std::string>>);
static_assert(container_ok<bounded_stack<int>>);
static_assert(container_ok<bounded_stack<std::string>>);
static_assert(wala::WALA_DEBUG || sizeof(vec<int>::iterator) == sizeof(int*));
static_assert(std::is_constructible_v<vec<int>, int, wala::uninit_t>);
static_assert(!std::is_constructible_v<vec<std::string>, int, wala::uninit_t>);
static_assert(!std::is_convertible_v<int, vec<int>>);
static_assert(!std::is_convertible_v<with_capacity, bounded_vec<int>>);

TEST_CASE("vec", "[wala::vec]") {
	vec<int> a(3);
	REQUIRE(a == vec<int>{0, 0, 0});
	vec<int> b(3, 7);
	REQUIRE(b == vec<int>{7, 7, 7});
	vec<int> c(4, uninit);
	REQUIRE(c.size() == 4);
	std::ranges::fill(c, 1);
	REQUIRE(c == vec<int>(4, 1));

	vec<int> d = {1, 2, 3};
	vec<int> e(std::from_range, std::vector{1, 2, 3});
	REQUIRE(d == e);
	REQUIRE(d != a);
	REQUIRE(d.front() == 1);
	REQUIRE(d.back() == 3);
	REQUIRE(std::as_const(d)[1] == 2);

	vec<int> f = d.clone();
	f[0] = 9;
	REQUIRE(d[0] == 1);
	REQUIRE(f != d);

	vec<int> g = std::move(f);
	REQUIRE(f.empty());
	REQUIRE(f.data() == nullptr);
	REQUIRE(f.begin() == f.end());
	REQUIRE(g == vec<int>{9, 2, 3});
	f = std::move(g);
	REQUIRE(f.size() == 3);

	std::span<int> sp(d);
	std::span<const int> csp(std::as_const(d));
	REQUIRE(sp.data() == d.data());
	REQUIRE(csp.size() == 3);
	std::ranges::sort(f, std::greater<>());
	REQUIRE(f == vec<int>{9, 3, 2});
	std::sort(f.begin(), f.end());
	REQUIRE(f == vec<int>{2, 3, 9});
	std::vector<int> sv(f.begin(), f.end());
	REQUIRE(sv == std::vector{2, 3, 9});
	REQUIRE(std::ranges::to<std::vector<int>>(f) == sv);
	REQUIRE(std::ranges::to<vec<int>>(std::views::iota(0, 4)) == vec<int>{0, 1, 2, 3});

	auto it = d.begin();
	it += 2;
	REQUIRE(*it == 3);
	REQUIRE(it[-1] == 2);
	REQUIRE(it - d.begin() == 2);
	REQUIRE(d.end() - it == 1);
	REQUIRE(--it == d.begin() + 1);
	REQUIRE(it++ == 1 + d.begin());
	REQUIRE(it < d.end());
	REQUIRE(std::to_address(d.end()) == d.data() + 3);
	vec<int>::const_iterator cit = it;
	REQUIRE(cit == std::as_const(d).begin() + 2);
	int sum = 0;
	for (int x : std::as_const(d)) sum += x;
	REQUIRE(sum == 6);

	vec<std::string> s(2, "ab");
	s.back() += "c";
	REQUIRE(s == vec<std::string>{"ab", "abc"});
	REQUIRE(s.clone() == s);

	bounded_vec<int> bv = vec<int>{1, 2, 3}.into_bounded();
	REQUIRE(bv.full());
	REQUIRE(bv.capacity() == 3);
	bounded_stack<int> bs = std::move(bv).into_stack();
	REQUIRE(bv.data() == nullptr);
	REQUIRE(bs.size() == 3);
	vec<int> v2 = std::move(bs).into_vec();
	REQUIRE(bs.data() == nullptr);
	REQUIRE(v2 == vec<int>{1, 2, 3});
	bs = std::move(v2).into_stack();
	REQUIRE(bs.full());
	REQUIRE(bs == bounded_stack<int>{1, 2, 3});
}

TEST_CASE("bounded_vec", "[wala::bounded_vec]") {
	bounded_vec<int> a(with_capacity{4});
	REQUIRE(a.empty());
	REQUIRE(a.capacity() == 4);
	a.push_back(1);
	REQUIRE(a.emplace_back(2) == 2);
	REQUIRE(a.size() == 2);
	REQUIRE(!a.full());
	a.grow_to(4, 9);
	REQUIRE(a.full());
	REQUIRE(a == bounded_vec<int>{1, 2, 9, 9});
	a.shrink_to(1);
	REQUIRE(a == bounded_vec<int>{1});
	a.pop_back();
	REQUIRE(a.empty());
	a.clear_and_set(3, 5);
	REQUIRE(a == bounded_vec<int>{5, 5, 5});
	a.grow_to(4);
	REQUIRE(a.back() == 0);
	a.clear();
	REQUIRE(a.empty());
	REQUIRE(a.capacity() == 4);
	a.grow_to(2, uninit);
	REQUIRE(a.size() == 2);

	bounded_vec<int> b(std::from_range, std::vector{1, 2, 3});
	REQUIRE(b.full());
	bounded_vec<int> c(std::from_range, std::vector{1, 2, 3}, with_capacity{5});
	REQUIRE(c.size() == 3);
	REQUIRE(c.capacity() == 5);
	REQUIRE(b == c);
	bounded_vec<int> d(std::from_range, std::views::iota(0, 10) | std::views::filter([](int x) { return x % 2 == 0; }), with_capacity{5});
	REQUIRE(d == bounded_vec<int>{0, 2, 4, 6, 8});
	auto e = std::ranges::to<bounded_vec<int>>(std::views::iota(0, 3), with_capacity{8});
	REQUIRE(e.capacity() == 8);
	REQUIRE(e == bounded_vec<int>{0, 1, 2});
	auto e2 = std::views::iota(0, 3) | std::ranges::to<bounded_vec<int>>();
	REQUIRE(e2.capacity() == 3);
	bounded_vec<int> f(2, 7, with_capacity{3});
	REQUIRE(f == bounded_vec<int>{7, 7});
	bounded_vec<int> g(2, uninit, with_capacity{3});
	REQUIRE(g.size() == 2);
	bounded_vec<int> h({1, 2}, with_capacity{3});
	REQUIRE(h.capacity() == 3);

	bounded_vec<int> cc = c.clone();
	REQUIRE(cc.capacity() == 5);
	cc[0] = 7;
	REQUIRE(c[0] == 1);
	bounded_vec<int> m = std::move(cc);
	REQUIRE(cc.capacity() == 0);
	REQUIRE(m[0] == 7);

	std::span<int> sp(c);
	REQUIRE(sp.size() == 3);
	std::ranges::sort(m, std::greater<>());
	REQUIRE(m == bounded_vec<int>{7, 3, 2});
	REQUIRE(m.end() - m.begin() == 3);
	REQUIRE(std::to_address(m.end()) == m.data() + 3);

	vec<int> v = std::move(b).into_vec();
	REQUIRE(v == vec<int>{1, 2, 3});
	REQUIRE(b.data() == nullptr);
	bounded_stack<int> s = std::move(c).into_stack();
	REQUIRE(s.size() == 3);
	REQUIRE(s.capacity() == 5);

	bounded_vec<std::string> str(with_capacity{2});
	str.emplace_back(3, 'x');
	str.push_back("y");
	REQUIRE(str.front() == "xxx");
	REQUIRE(str.clone() == str);}

TEST_CASE("bounded_stack", "[wala::bounded_stack]") {
	bounded_stack<int> a(with_capacity{4});
	REQUIRE(a.empty());
	REQUIRE(a.capacity() == 4);
	a.push_back(1);
	REQUIRE(a.emplace_back(2) == 2);
	REQUIRE(a.size() == 2);
	REQUIRE(!a.full());
	REQUIRE(std::to_address(a.end()) == a.top);
	a.grow_to(4, 9);
	REQUIRE(a.full());
	REQUIRE(a == bounded_stack<int>{1, 2, 9, 9});
	a.shrink_to(1);
	REQUIRE(a == bounded_stack<int>{1});
	a.pop_back();
	REQUIRE(a.empty());
	a.clear_and_set(3, 5);
	REQUIRE(a == bounded_stack<int>{5, 5, 5});
	a.grow_to(4);
	REQUIRE(a.back() == 0);
	a.clear();
	REQUIRE(a.empty());
	a.grow_to(2, uninit);
	REQUIRE(a.size() == 2);

	bounded_stack<int> b(std::from_range, std::vector{1, 2, 3});
	REQUIRE(b.full());
	bounded_stack<int> c(std::from_range, std::vector{1, 2, 3}, with_capacity{5});
	REQUIRE(c.size() == 3);
	REQUIRE(c.capacity() == 5);
	REQUIRE(b == c);
	bounded_stack<int> d(std::from_range, std::views::iota(0, 10) | std::views::filter([](int x) { return x % 2 == 0; }), with_capacity{5});
	REQUIRE(d == bounded_stack<int>{0, 2, 4, 6, 8});
	auto e = std::ranges::to<bounded_stack<int>>(std::views::iota(0, 3), with_capacity{8});
	REQUIRE(e.capacity() == 8);
	bounded_stack<int> f(2, 7, with_capacity{3});
	REQUIRE(f == bounded_stack<int>{7, 7});
	bounded_stack<int> g(2, uninit, with_capacity{3});
	REQUIRE(g.size() == 2);
	bounded_stack<int> h({1, 2}, with_capacity{3});
	REQUIRE(h.capacity() == 3);

	bounded_stack<int> cc = c.clone();
	REQUIRE(cc.capacity() == 5);
	cc[0] = 7;
	REQUIRE(c[0] == 1);
	bounded_stack<int> m = std::move(cc);
	REQUIRE(cc.capacity() == 0);
	REQUIRE(m[0] == 7);

	std::span<int> sp(c);
	REQUIRE(sp.size() == 3);
	std::ranges::sort(m, std::greater<>());
	REQUIRE(m == bounded_stack<int>{7, 3, 2});
	REQUIRE(m.end() - m.begin() == 3);

	vec<int> v = std::move(b).into_vec();
	REQUIRE(v == vec<int>{1, 2, 3});
	REQUIRE(b.data() == nullptr);
	bounded_vec<int> bv = std::move(c).into_bounded();
	REQUIRE(bv.size() == 3);
	REQUIRE(bv.capacity() == 5);
	REQUIRE(c.data() == nullptr);

	bounded_stack<std::string> str(with_capacity{2});
	str.emplace_back(3, 'x');
	str.push_back("y");
	REQUIRE(str.front() == "xxx");
	REQUIRE(str.clone() == str);
}
