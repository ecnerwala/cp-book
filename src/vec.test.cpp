#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include "vec.hpp"

#include <array>
#include <memory>
#include <ranges>
#include <string>
#include <vector>

#include <csignal>
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

template <class C> concept contiguous_container =
	std::ranges::contiguous_range<C>
	&& std::ranges::contiguous_range<const C>
	&& std::ranges::sized_range<C>
	&& std::contiguous_iterator<typename C::iterator>
	&& std::contiguous_iterator<typename C::const_iterator>
	&& !std::convertible_to<typename C::const_iterator, typename C::iterator>
	&& std::convertible_to<typename C::iterator, typename C::const_iterator>
	&& std::same_as<std::iter_value_t<typename C::const_iterator>, typename C::value_type>
	&& std::same_as<std::ranges::iterator_t<C>, typename C::iterator>
	&& std::same_as<std::ranges::iterator_t<const C>, typename C::const_iterator>
	&& std::same_as<std::ranges::range_reference_t<C>, typename C::value_type&>
	&& std::same_as<std::ranges::range_reference_t<const C>, const typename C::value_type&>
	&& std::same_as<std::ranges::range_size_t<C>, int>
	&& std::same_as<std::ranges::range_difference_t<C>, int>
	&& std::is_nothrow_move_constructible_v<C> && std::is_nothrow_move_assignable_v<C> && std::is_nothrow_swappable_v<C>
	&& std::copyable<C> == std::copy_constructible<typename C::value_type>
	&& std::equality_comparable<C>;

static_assert(contiguous_container<wala::vec<bool>>);
static_assert(contiguous_container<wala::vec<int>>);
static_assert(contiguous_container<wala::vec<std::string>>);
static_assert(contiguous_container<wala::vec<std::unique_ptr<std::string>>>);
static_assert(contiguous_container<wala::bounded_vec<bool>>);
static_assert(contiguous_container<wala::bounded_vec<int>>);
static_assert(contiguous_container<wala::bounded_vec<std::string>>);
static_assert(contiguous_container<wala::bounded_vec<std::unique_ptr<std::string>>>);
static_assert(contiguous_container<wala::bounded_stack<bool>>);
static_assert(contiguous_container<wala::bounded_stack<int>>);
static_assert(contiguous_container<wala::bounded_stack<std::string>>);
static_assert(contiguous_container<wala::bounded_stack<std::unique_ptr<std::string>>>);
static_assert(contiguous_container<wala::growable_vec<bool>>);
static_assert(contiguous_container<wala::growable_vec<int>>);
static_assert(contiguous_container<wala::growable_vec<std::string>>);
static_assert(contiguous_container<wala::growable_vec<std::unique_ptr<std::string>>>);

TEST_CASE("braced fill values", "[wala::vec]") {
	using P = std::array<int, 2>;
	wala::vec<P> a(2, {-1, -1});
	REQUIRE(a == wala::vec<P>{{-1, -1}, {-1, -1}});
	wala::vec<P> z(2, {});
	REQUIRE(z == wala::vec<P>{{0, 0}, {0, 0}});
	wala::vec<std::string> s(2, "x");
	REQUIRE(s == wala::vec<std::string>{"x", "x"});

	wala::bounded_vec<P> b(wala::with_capacity{6}, 2, {1, 2});
	REQUIRE(b == wala::bounded_vec<P>{{1, 2}, {1, 2}});
	REQUIRE(b.capacity() == 6);
	b.grow_by(1, {3, 4});
	b.grow_to(4, {});
	b.insert_n(0, 1, {5, 6});
	REQUIRE(b == wala::bounded_vec<P>{{5, 6}, {1, 2}, {1, 2}, {3, 4}, {0, 0}});
	b.grow_to(b.end() + 1, {7, 8});
	REQUIRE(b.back() == P{7, 8});
	wala::bounded_vec<P> b2(3, {9, 9});
	REQUIRE(b2.capacity() == 3);
	REQUIRE(b2[2] == P{9, 9});

	wala::bounded_stack<P> c(wala::with_capacity{4}, 1, {1, 1});
	c.grow_by(1, {2, 2});
	c.grow_to(3, {3, 3});
	c.grow_to(c.end() + 1, {4, 4});
	REQUIRE(c == wala::bounded_stack<P>{{1, 1}, {2, 2}, {3, 3}, {4, 4}});
	wala::bounded_stack<P> c2(2, {});
	REQUIRE(c2 == wala::bounded_stack<P>{{0, 0}, {0, 0}});

	wala::growable_vec<P> g(1, {1, 1});
	g.grow_by(2, {2, 2});
	g.grow_to(5, {});
	g.insert_n(0, 1, {3, 3});
	g.push_front({4, 4});
	REQUIRE(g == wala::growable_vec<P>{{4, 4}, {3, 3}, {1, 1}, {2, 2}, {2, 2}, {0, 0}, {0, 0}});
}

// Only the kinds of construction that make sense for T are accepted
static_assert(std::constructible_from<wala::vec<int>, int, wala::uninit_t>);
static_assert(!std::constructible_from<wala::vec<std::string>, int, wala::uninit_t>);
static_assert(std::constructible_from<wala::vec<std::string>, int, const char*>);
static_assert(!std::constructible_from<wala::vec<std::unique_ptr<int>>, std::initializer_list<std::unique_ptr<int>>>);
static_assert(std::constructible_from<wala::bounded_vec<int>, wala::with_capacity, int, wala::uninit_t>);
static_assert(!std::constructible_from<wala::bounded_stack<std::string>, wala::with_capacity, int, wala::uninit_t>);

namespace {

// Tracks live instances so tests can check construction/destruction balance.
struct tracked {
	static inline int live = 0;
	std::string s;
	tracked() : s("") { live++; }
	tracked(std::string s_) : s(std::move(s_)) { live++; }
	tracked(const tracked& o) : s(o.s) { live++; }
	tracked(tracked&& o) noexcept : s(std::move(o.s)) { live++; }
	tracked& operator=(const tracked&) = default;
	tracked& operator=(tracked&&) = default;
	~tracked() { live--; }
	friend bool operator==(const tracked& a, const tracked& b) { return a.s == b.s; }
};

template <class C> std::string show(const C& c) {
	std::string r;
	for (const auto& x : c) {
		r += x.s;
		r += ',';
	}
	return r;
}

// Runs f in a forked child and checks that it died from a failed debug_assert (abort).
template <class F> void require_debug_abort(F f) {
	if constexpr (wala::WALA_DEBUG) {
		pid_t pid = fork();
		REQUIRE(pid >= 0);
		if (pid == 0) {
			// Keep Catch2's signal handlers and reporters out of the child
			signal(SIGABRT, SIG_DFL);
			int devnull = open("/dev/null", O_WRONLY);
			dup2(devnull, STDOUT_FILENO);
			dup2(devnull, STDERR_FILENO);
			f();
			_exit(0);
		}
		int status = 0;
		REQUIRE(waitpid(pid, &status, 0) == pid);
		REQUIRE(WIFSIGNALED(status));
		REQUIRE(WTERMSIG(status) == SIGABRT);
	}
}

} // namespace

TEST_CASE("vec construction", "[wala::vec]") {
	{
		wala::vec<int> a;
		REQUIRE(a.empty());
		REQUIRE(a.size() == 0);
		REQUIRE(a.begin() == a.end());
		REQUIRE(a.data() == nullptr);

		wala::vec<int> b(3);
		REQUIRE(b == wala::vec<int>{0, 0, 0});
		wala::vec<int> c(3, 7);
		REQUIRE(c == wala::vec<int>{7, 7, 7});
		wala::vec<int> d(3, wala::uninit);
		REQUIRE(d.size() == 3);
		for (int i = 0; i < d.size(); i++) d[i] = i;
		REQUIRE(d == wala::vec<int>{0, 1, 2});

		std::vector<int> src{1, 2, 3, 4};
		wala::vec<int> e(std::from_range, src);
		REQUIRE(std::ranges::equal(e, src));
		wala::vec<int> f(std::from_range, src | std::views::transform([](int x) { return x * x; }));
		REQUIRE(f == wala::vec<int>{1, 4, 9, 16});
		wala::vec<int> g(std::from_range, std::views::iota(0, 0));
		REQUIRE(g.empty());
	}
	{
		wala::vec<tracked> a(2, tracked("x"));
		REQUIRE(tracked::live == 2);
		REQUIRE(show(a) == "x,x,");
		wala::vec<tracked> b(3);
		REQUIRE(tracked::live == 5);
		REQUIRE(show(b) == ",,,");
		wala::vec<tracked> c(std::from_range, std::vector<std::string>{"p", "q"});
		REQUIRE(tracked::live == 7);
		REQUIRE(show(c) == "p,q,");
	}
	REQUIRE(tracked::live == 0);
}

TEST_CASE("vec access and comparison", "[wala::vec]") {
	wala::vec<std::string> a{"a", "b", "c"};
	REQUIRE(a.size() == 3);
	REQUIRE(a.front() == "a");
	REQUIRE(a.back() == "c");
	REQUIRE(a[1] == "b");
	REQUIRE(a.data()[2] == "c");
	REQUIRE(*a.begin() == "a");
	REQUIRE(*(a.end() - 1) == "c");
	REQUIRE(a.end() - a.begin() == 3);
	REQUIRE(a.begin()[1] == "b");
	REQUIRE(std::ranges::distance(a) == 3);
	a[1] = "B";
	REQUIRE(a == wala::vec<std::string>{"a", "B", "c"});
	REQUIRE(a != wala::vec<std::string>{"a", "b", "c"});
	REQUIRE(a < wala::vec<std::string>{"a", "C"});
	REQUIRE(a > wala::vec<std::string>{"a"});
	REQUIRE((a <=> wala::vec<std::string>{"a", "B", "c"}) == 0);

	const wala::vec<std::string>& ca = a;
	REQUIRE(ca[0] == "a");
	REQUIRE(ca.data() == a.data());
	REQUIRE(std::to_address(a.begin()) == a.data());
	REQUIRE(std::to_address(ca.end()) == ca.data() + 3);
	REQUIRE(ca.front() == "a");
	REQUIRE(ca.back() == "c");
	REQUIRE(*ca.begin() == "a");
	static_assert(std::same_as<decltype(ca.begin()), wala::vec<std::string>::const_iterator>);
	wala::vec<std::string>::const_iterator it = a.begin();
	REQUIRE(it == ca.begin());
	REQUIRE(a.begin() == ca.begin());

	// Standard algorithms and conversions through the iterators
	std::vector<std::string> v(a.begin(), a.end());
	REQUIRE(v == std::vector<std::string>{"a", "B", "c"});
	std::span<std::string> sp(a);
	REQUIRE(sp.size() == 3);
	REQUIRE(std::ranges::to<std::vector<std::string>>(a) == v);
	std::ranges::sort(a);
	REQUIRE(a == wala::vec<std::string>{"B", "a", "c"});
}

TEST_CASE("vec copy and move", "[wala::vec]") {
	{
		wala::vec<tracked> a(std::from_range, std::vector<std::string>{"a", "b"});
		wala::vec<tracked> b = a;
		REQUIRE(b == a);
		REQUIRE(b.data() != a.data());
		REQUIRE(tracked::live == 4);
		b[0].s = "z";
		REQUIRE(a[0].s == "a");

		wala::vec<tracked> c = std::move(a);
		REQUIRE(show(c) == "a,b,");
		REQUIRE(a.empty());
		REQUIRE(tracked::live == 4);

		a = c;
		REQUIRE(a == c);
		REQUIRE(tracked::live == 6);
		a = std::move(b);
		REQUIRE(show(a) == "z,b,");
		// Moved-from is valid but unspecified
		b = wala::vec<tracked>();
		REQUIRE(b.empty());
		REQUIRE(tracked::live == 4);
		auto& self = a;
		a = self;
		REQUIRE(show(a) == "z,b,");
		swap(a, c);
		REQUIRE(show(a) == "a,b,");
		REQUIRE(show(c) == "z,b,");
	}
	REQUIRE(tracked::live == 0);

	wala::vec<wala::vec<int>> nested{{1, 2}, {3}};
	wala::vec<wala::vec<int>> nested2 = nested;
	nested2[0][0] = 9;
	REQUIRE(nested[0][0] == 1);
	REQUIRE(nested2 == wala::vec<wala::vec<int>>{{9, 2}, {3}});

	wala::vec<std::unique_ptr<int>> ptrs(2);
	ptrs[0] = std::make_unique<int>(5);
	wala::vec<std::unique_ptr<int>> ptrs2 = std::move(ptrs);
	REQUIRE(*ptrs2[0] == 5);
	REQUIRE(ptrs2[1] == nullptr);
	wala::vec<std::unique_ptr<int>> ptrs3(std::from_range, std::views::as_rvalue(ptrs2));
	REQUIRE(*ptrs3[0] == 5);
	REQUIRE(ptrs2[0] == nullptr);
}

TEST_CASE("vec debug checks", "[wala::vec]") {
	wala::vec<int> a{1, 2, 3};
	require_debug_abort([&] { (void)a[3]; });
	require_debug_abort([&] { (void)a[-1]; });
	require_debug_abort([&] { (void)*a.end(); });
	require_debug_abort([&] { (void)a.begin()[3]; });
	require_debug_abort([&] { wala::vec<int> b; (void)b.front(); });
	require_debug_abort([&] { wala::vec<int> b; (void)b.back(); });
	require_debug_abort([&] { wala::vec<int> b(-1); });
}

template <class C, class... A> C make_bounded(int cap, A&&... a) {
	return C(wala::with_capacity{cap}, std::forward<A>(a)...);
}

TEMPLATE_TEST_CASE("bounded construction", "[wala::vec]", wala::bounded_vec<int>, wala::bounded_stack<int>, wala::growable_vec<int>) {
	using C = TestType;
	C a;
	REQUIRE(a.empty());
	REQUIRE(a.size() == 0);
	REQUIRE(a.capacity() == 0);
	REQUIRE(a.data() == nullptr);

	C b(wala::with_capacity{5});
	REQUIRE(b.empty());
	REQUIRE(b.capacity() == 5);
	REQUIRE(b.data() != nullptr);

	C c(wala::with_capacity{5}, 3);
	REQUIRE(c.capacity() == 5);
	REQUIRE(c == C{0, 0, 0});
	C d(wala::with_capacity{5}, 3, 7);
	REQUIRE(d == C{7, 7, 7});
	C e(wala::with_capacity{5}, 3, wala::uninit);
	REQUIRE(e.size() == 3);
	REQUIRE(e.capacity() == 5);

	// Without with_capacity, capacity == size
	C f(3, 1);
	REQUIRE(f.capacity() == 3);
	REQUIRE(f == C{1, 1, 1});
	C g{4, 5};
	REQUIRE(g.capacity() == 2);
	REQUIRE(g.size() == 2);
	REQUIRE(g[1] == 5);
	C h(wala::with_capacity{4}, {4, 5});
	REQUIRE(h.capacity() == 4);
	REQUIRE(h == g);

	std::vector<int> src{1, 2, 3};
	C i(std::from_range, src);
	REQUIRE(i.capacity() == 3);
	REQUIRE(std::ranges::equal(i, src));
	C j(wala::with_capacity{6}, std::from_range, src);
	REQUIRE(j.capacity() == 6);
	REQUIRE(j == i);
	// Unsized input ranges only need the explicit capacity
	C k(wala::with_capacity{6}, std::from_range, src | std::views::filter([](int x) { return x != 2; }));
	REQUIRE(k == C{1, 3});
	REQUIRE(k.capacity() == 6);
}

TEMPLATE_TEST_CASE("bounded push/pop", "[wala::vec]", wala::bounded_vec<tracked>, wala::bounded_stack<tracked>, wala::growable_vec<tracked>) {
	using C = TestType;
	{
		C a(wala::with_capacity{4});
		REQUIRE(a.empty());
		tracked& r = a.emplace_back("a");
		REQUIRE(r.s == "a");
		REQUIRE(&r == &a.back());
		a.push_back(tracked("b"));
		tracked c("c");
		a.push_back(c);
		REQUIRE(tracked::live == 4);
		REQUIRE(show(a) == "a,b,c,");
		REQUIRE(a.size() == 3);
		REQUIRE(a.capacity() == 4);
		// Pushing an element of the container itself
		a.push_back(a.front());
		REQUIRE(show(a) == "a,b,c,a,");
		a.pop_back();
		REQUIRE(show(a) == "a,b,c,");
		REQUIRE(tracked::live == 4);
		a.pop_back();
		a.pop_back();
		a.pop_back();
		REQUIRE(a.empty());
		REQUIRE(a.capacity() == 4);
		REQUIRE(tracked::live == 1);
		// Reusing the slots after popping
		a.emplace_back("d");
		REQUIRE(show(a) == "d,");
		// Keep c alive to here
		REQUIRE(c.s == "c");
	}
	REQUIRE(tracked::live == 0);
}

TEMPLATE_TEST_CASE("bounded grow/shrink", "[wala::vec]", wala::bounded_vec<tracked>, wala::bounded_stack<tracked>, wala::growable_vec<tracked>) {
	using C = TestType;
	{
		C a(wala::with_capacity{10});
		auto [f1, l1] = a.grow_by(2, tracked("x"));
		REQUIRE(f1 == a.begin());
		REQUIRE(l1 == a.end());
		REQUIRE(l1 - f1 == 2);
		REQUIRE(show(a) == "x,x,");

		auto [f2, l2] = a.grow_to(5);
		REQUIRE(f2 == a.begin() + 2);
		REQUIRE(l2 == a.end());
		REQUIRE(show(a) == "x,x,,,,");
		REQUIRE(tracked::live == 5);

		// Growing by a value from the container itself
		auto [f3, l3] = a.grow_by(1, a.front());
		REQUIRE(f3->s == "x");
		REQUIRE(show(a) == "x,x,,,,x,");

		auto r0 = a.grow_by(0);
		REQUIRE(r0.begin() == a.end());
		REQUIRE(r0.empty());
		auto r1 = a.grow_to(6, tracked("unused"));
		REQUIRE(r1.empty());

		a.shrink_by(1);
		REQUIRE(show(a) == "x,x,,,,");
		a.shrink_to(2);
		REQUIRE(show(a) == "x,x,");
		REQUIRE(tracked::live == 2);
		a.shrink_by(0);
		REQUIRE(a.size() == 2);
		REQUIRE(a.capacity() == 10);
		a.clear();
		REQUIRE(a.empty());
		REQUIRE(a.capacity() == 10);
		REQUIRE(tracked::live == 0);
	}
	REQUIRE(tracked::live == 0);
	{
		wala::bounded_vec<int> b(wala::with_capacity{4}, 1, 3);
		auto [f, l] = b.grow_by(2, wala::uninit);
		REQUIRE(l - f == 2);
		f[0] = 4;
		f[1] = 5;
		REQUIRE(b == wala::bounded_vec<int>{3, 4, 5});
	}
}

TEMPLATE_TEST_CASE("bounded iterator marks", "[wala::vec]", wala::bounded_vec<tracked>, wala::bounded_stack<tracked>) {
	using C = TestType;
	{
		C a(wala::with_capacity{10});
		a.emplace_back("a");
		auto mark = a.end();
		a.emplace_back("b");
		a.emplace_back("c");
		REQUIRE(tracked::live == 3);
		a.shrink_to(mark);
		REQUIRE(show(a) == "a,");
		REQUIRE(tracked::live == 1);
		REQUIRE(mark == a.end());
		auto [f, l] = a.grow_to(a.begin() + 4, tracked("z"));
		REQUIRE(f == mark);
		REQUIRE(l == a.end());
		REQUIRE(show(a) == "a,z,z,z,");
		a.shrink_to(a.begin());
		REQUIRE(a.empty());
		a.grow_to(a.end());
		REQUIRE(a.empty());
	}
	REQUIRE(tracked::live == 0);
}

TEMPLATE_TEST_CASE("bounded append_range", "[wala::vec]", wala::bounded_vec<tracked>, wala::bounded_stack<tracked>, wala::growable_vec<tracked>) {
	using C = TestType;
	{
		C a(wala::with_capacity{8});
		auto [f1, l1] = a.append_range(std::vector<std::string>{"a", "b"});
		REQUIRE(f1 == a.begin());
		REQUIRE(l1 == a.end());
		REQUIRE(show(a) == "a,b,");
		std::vector<tracked> src{tracked("c"), tracked("d")};
		auto [f2, l2] = a.append_range(src);
		REQUIRE(f2 == a.begin() + 2);
		REQUIRE(l2 == a.end());
		REQUIRE(show(a) == "a,b,c,d,");
		REQUIRE(src[0].s == "c");
		// Moving out of an rvalue range
		a.append_range(src | std::views::as_rvalue);
		REQUIRE(show(a) == "a,b,c,d,c,d,");
		REQUIRE(src[0].s == "");
		REQUIRE(src[1].s == "");
		// Unsized input range, and an empty one
		auto r = a.append_range(std::views::iota(0, 4) | std::views::filter([](int x) { return x % 2; }) | std::views::transform([](int x) { return std::to_string(x); }));
		REQUIRE(r.size() == 2);
		REQUIRE(show(a) == "a,b,c,d,c,d,1,3,");
		REQUIRE(a.append_range(std::vector<tracked>{}).empty());
		REQUIRE(a.size() == 8);
		REQUIRE(tracked::live == 10);
	}
	REQUIRE(tracked::live == 0);
}

TEMPLATE_TEST_CASE("bounded copy and move", "[wala::vec]", wala::bounded_vec<tracked>, wala::bounded_stack<tracked>, wala::growable_vec<tracked>) {
	using C = TestType;
	{
		C a(wala::with_capacity{5}, std::from_range, std::vector<std::string>{"a", "b"});
		C b = a;
		REQUIRE(b == a);
		REQUIRE(b.capacity() == 5);
		REQUIRE(b.data() != a.data());
		REQUIRE(tracked::live == 4);

		C c = std::move(a);
		REQUIRE(show(c) == "a,b,");
		REQUIRE(c.capacity() == 5);
		REQUIRE(a.empty());
		REQUIRE(a.capacity() == 0);
		REQUIRE(tracked::live == 4);

		a = c;
		REQUIRE(a == c);
		REQUIRE(a.capacity() == 5);
		a = std::move(b);
		REQUIRE(show(a) == "a,b,");
		// Moved-from is valid but unspecified
		b = C();
		REQUIRE(b.empty());
		REQUIRE(tracked::live == 4);
		auto& self = a;
		a = self;
		REQUIRE(show(a) == "a,b,");
		swap(a, c);
		REQUIRE(a == c);
	}
	REQUIRE(tracked::live == 0);

	C ptrs(wala::with_capacity{2});
	static_assert(!std::copyable<wala::bounded_vec<std::unique_ptr<int>>>);
	static_assert(!std::copyable<wala::bounded_stack<std::unique_ptr<int>>>);
	wala::bounded_stack<std::unique_ptr<int>> p(wala::with_capacity{2});
	p.emplace_back(std::make_unique<int>(3));
	p.push_back(std::make_unique<int>(4));
	wala::bounded_stack<std::unique_ptr<int>> q = std::move(p);
	REQUIRE(*q[0] == 3);
	REQUIRE(*q.back() == 4);
}

TEMPLATE_TEST_CASE("bounded front ops", "[wala::vec]", wala::bounded_vec<tracked>, wala::bounded_stack<tracked>, wala::growable_vec<tracked>) {
	using C = TestType;
	{
		C a(wala::with_capacity{6}, std::from_range, std::vector<std::string>{"b", "c"});
		a.emplace_front("a");
		REQUIRE(show(a) == "a,b,c,");
		a.push_front(tracked("z"));
		REQUIRE(show(a) == "z,a,b,c,");
		// Pushing an element of the container itself
		a.push_front(a.back());
		REQUIRE(show(a) == "c,z,a,b,c,");
		REQUIRE(tracked::live == 5);
		a.pop_front();
		a.pop_front();
		REQUIRE(show(a) == "a,b,c,");
		REQUIRE(tracked::live == 3);
	}
	REQUIRE(tracked::live == 0);
}

TEMPLATE_TEST_CASE("bounded insert_n/erase_n", "[wala::vec]", wala::bounded_vec<tracked>, wala::bounded_stack<tracked>, wala::growable_vec<tracked>) {
	using C = TestType;
	{
		C a(wala::with_capacity{8}, std::from_range, std::vector<std::string>{"a", "b"});
		a.insert_n(0, 2, tracked("0"));
		REQUIRE(show(a) == "0,0,a,b,");
		a.insert_n(3, 1, a.back());
		REQUIRE(show(a) == "0,0,a,b,b,");
		a.insert_n(5, 2);
		REQUIRE(show(a) == "0,0,a,b,b,,,");
		REQUIRE(tracked::live == 7);
		a.erase_n(0, 2);
		REQUIRE(show(a) == "a,b,b,,,");
		a.erase_n(2, 3);
		REQUIRE(show(a) == "a,b,");
		a.erase_n(2, 0);
		a.erase_n(0, 0);
		REQUIRE(show(a) == "a,b,");
		REQUIRE(tracked::live == 2);
	}
	REQUIRE(tracked::live == 0);
}

TEMPLATE_TEST_CASE("bounded debug checks", "[wala::vec]", wala::bounded_vec<int>, wala::bounded_stack<int>) {
	using C = TestType;
	C a(wala::with_capacity{3}, {1, 2});
	require_debug_abort([&] { (void)a[2]; });
	require_debug_abort([&] { (void)*a.end(); });
	require_debug_abort([&] { a.emplace_back(3); a.emplace_back(4); });
	require_debug_abort([&] { a.grow_by(2); });
	require_debug_abort([&] { a.grow_by(-1); });
	require_debug_abort([&] { a.grow_to(1); });
	require_debug_abort([&] { a.grow_to(4); });
	require_debug_abort([&] { a.shrink_by(3); });
	require_debug_abort([&] { a.shrink_to(3); });
	require_debug_abort([&] { a.shrink_to(-1); });
	require_debug_abort([&] { a.append_range(std::vector<int>{3, 4}); });
	require_debug_abort([&] { C b(wala::with_capacity{2}, 3); });
	require_debug_abort([&] { C b(wala::with_capacity{1}, std::from_range, std::vector<int>{1, 2}); });
	require_debug_abort([&] { C b; (void)b.back(); });
	require_debug_abort([&] { C b(wala::with_capacity{1}); b.pop_back(); });
	require_debug_abort([&] { a.erase_n(1, 2); });
	require_debug_abort([&] { a.insert_n(3, 1); });
}

TEST_CASE("growable_vec growth", "[wala::vec]") {
	{
		wala::growable_vec<tracked> a;
		REQUIRE(a.capacity() == 0);
		int reallocs = 0;
		for (int i = 0; i < 100; i++) {
			const tracked* old_data = a.data();
			int old_cap = a.capacity();
			tracked& r = a.emplace_back(std::to_string(i));
			REQUIRE(&r == &a.back());
			REQUIRE(a.size() == i + 1);
			REQUIRE(a.capacity() >= a.size());
			if (a.data() != old_data) {
				reallocs++;
				// Doubling
				REQUIRE(a.capacity() == std::max(2 * old_cap, 1));
			}
		}
		REQUIRE(reallocs == 8);
		REQUIRE(a.capacity() == 128);
		REQUIRE(tracked::live == 100);
		for (int i = 0; i < 100; i++) REQUIRE(a[i].s == std::to_string(i));

		// Returned blocks are valid after reallocation
		auto [f, l] = a.grow_by(100, tracked("x"));
		REQUIRE(a.capacity() == 256);
		REQUIRE(f == a.begin() + 100);
		REQUIRE(l == a.end());
		REQUIRE(f->s == "x");
		REQUIRE(a[99].s == "99");
		auto [f2, l2] = a.grow_to(1000);
		REQUIRE(a.capacity() == 1000);
		REQUIRE(f2 == a.begin() + 200);
		REQUIRE(l2 == a.end());
		REQUIRE(tracked::live == 1000);

		// Growing by an element of the container itself across a reallocation
		a.shrink_to(1000);
		a.shrink_to_fit();
		REQUIRE(a.capacity() == 1000);
		a.push_back(a.front());
		REQUIRE(a.capacity() == 2000);
		REQUIRE(a.back().s == "0");
		a.grow_by(1500, a[99]);
		REQUIRE(a.capacity() == 4000);
		REQUIRE(a.back().s == "99");
		REQUIRE(a.size() == 2501);

		a.shrink_to(3);
		REQUIRE(a.capacity() == 4000);
		a.shrink_to_fit();
		REQUIRE(a.capacity() == 3);
		REQUIRE(show(a) == "0,1,2,");
		a.set_capacity(10);
		REQUIRE(a.capacity() == 10);
		REQUIRE(show(a) == "0,1,2,");
		a.clear();
		a.shrink_to_fit();
		REQUIRE(a.capacity() == 0);
		REQUIRE(a.data() == nullptr);
		REQUIRE(tracked::live == 0);
	}
	REQUIRE(tracked::live == 0);
	{
		// Inherited ops grow through the overriding growers
		wala::growable_vec<tracked> a{tracked("b")};
		REQUIRE(a.capacity() == 1);
		a.push_front(tracked("a"));
		REQUIRE(show(a) == "a,b,");
		a.insert_n(2, 3, tracked("c"));
		REQUIRE(show(a) == "a,b,c,c,c,");
		a.insert(1, a.back());
		REQUIRE(show(a) == "a,c,b,c,c,c,");
		a.emplace_front("z");
		REQUIRE(show(a) == "z,a,c,b,c,c,c,");
		// Sized and unsized ranges
		a.append_range(std::vector<std::string>{"d", "e"});
		auto r = a.append_range(std::views::iota(0, 10) | std::views::filter([](int x) { return x % 5 == 0; }) | std::views::transform([](int x) { return std::to_string(x); }));
		REQUIRE(r.size() == 2);
		REQUIRE(r.begin() == a.end() - 2);
		REQUIRE(show(a) == "z,a,c,b,c,c,c,d,e,0,5,");
		REQUIRE(a.capacity() >= a.size());
		REQUIRE(tracked::live == 11);
		a.erase_n(0, 9);
		a.pop_front();
		REQUIRE(show(a) == "5,");
	}
	REQUIRE(tracked::live == 0);
	{
		// Capacity is preserved by copy and move, like bounded_vec
		wala::growable_vec<tracked> a(wala::with_capacity{8}, 2, tracked("a"));
		wala::growable_vec<tracked> b = a;
		REQUIRE(b.capacity() == 8);
		REQUIRE(b == a);
		wala::growable_vec<tracked> c = std::move(a);
		REQUIRE(c.capacity() == 8);
		REQUIRE(a.capacity() == 0);
		swap(b, c);
		REQUIRE(b == c);
		REQUIRE(tracked::live == 4);
	}
	REQUIRE(tracked::live == 0);

	wala::growable_vec<std::unique_ptr<int>> p;
	for (int i = 0; i < 10; i++) p.emplace_back(std::make_unique<int>(i));
	p.shrink_to_fit();
	REQUIRE(p.capacity() == 10);
	REQUIRE(*p[7] == 7);
	wala::growable_vec<wala::growable_vec<int>> nested{{1, 2}, {3}};
	nested.emplace_back(std::from_range, std::views::iota(4, 7));
	nested[0].push_back(9);
	wala::growable_vec<wala::growable_vec<int>> nested2 = nested;
	REQUIRE(nested2 == wala::growable_vec<wala::growable_vec<int>>{{1, 2, 9}, {3}, {4, 5, 6}});
}

TEST_CASE("growable_vec debug checks", "[wala::vec]") {
	wala::growable_vec<int> a{1, 2};
	require_debug_abort([&] { (void)a[2]; });
	require_debug_abort([&] { (void)*a.end(); });
	require_debug_abort([&] { a.grow_by(-1); });
	require_debug_abort([&] { a.grow_to(1); });
	require_debug_abort([&] { a.shrink_by(3); });
	require_debug_abort([&] { a.set_capacity(1); });
	require_debug_abort([&] { wala::growable_vec<int> b; (void)b.back(); });
	require_debug_abort([&] { wala::growable_vec<int> b; b.pop_back(); });
}

TEST_CASE("into_* conversions", "[wala::vec]") {
	using V = wala::vec<tracked>;
	using B = wala::bounded_vec<tracked>;
	using S = wala::bounded_stack<tracked>;
	using G = wala::growable_vec<tracked>;
	{
		V v(std::from_range, std::vector<std::string>{"a", "b", "c"});
		const tracked* p = v.data();
		B b = std::move(v).into_bounded();
		REQUIRE(v.empty());
		REQUIRE(v.data() == nullptr);
		REQUIRE(b.data() == p);
		REQUIRE(b.size() == 3);
		REQUIRE(b.capacity() == 3);
		REQUIRE(show(b) == "a,b,c,");
		REQUIRE(tracked::live == 3);

		S s = std::move(b).into_stack();
		REQUIRE(b.empty());
		REQUIRE(b.capacity() == 0);
		REQUIRE(b.data() == nullptr);
		REQUIRE(s.data() == p);
		REQUIRE(s.size() == 3);
		REQUIRE(s.capacity() == 3);
		REQUIRE(show(s) == "a,b,c,");

		G g = std::move(s).into_growable();
		REQUIRE(s.empty());
		REQUIRE(s.capacity() == 0);
		REQUIRE(s.data() == nullptr);
		REQUIRE(g.data() == p);
		REQUIRE(g.size() == 3);
		REQUIRE(g.capacity() == 3);
		g.emplace_back("d");
		REQUIRE(g.capacity() == 6);
		REQUIRE(show(g) == "a,b,c,d,");
		REQUIRE(tracked::live == 4);

		// Not full, so into_vec shrinks first
		const tracked* q = g.data();
		V v2 = std::move(g).into_vec();
		REQUIRE(g.empty());
		REQUIRE(g.capacity() == 0);
		REQUIRE(v2.data() != q);
		REQUIRE(show(v2) == "a,b,c,d,");
		REQUIRE(tracked::live == 4);

		g = std::move(v2).into_growable();
		REQUIRE(g.capacity() == 4);
		B b2 = std::move(g).into_bounded();
		REQUIRE(g.empty());
		REQUIRE(b2.capacity() == 4);
		G g2 = std::move(b2).into_growable();
		REQUIRE(b2.empty());
		REQUIRE(show(g2) == "a,b,c,d,");
		S s2 = std::move(g2).into_stack();
		REQUIRE(g2.empty());
		REQUIRE(s2.capacity() == 4);
		V v3 = std::move(s2).into_vec();
		REQUIRE(s2.empty());
		REQUIRE(v3.size() == 4);
		REQUIRE(show(v3) == "a,b,c,d,");
		REQUIRE(tracked::live == 4);
	}
	REQUIRE(tracked::live == 0);
	{
		// Full containers convert to vec in place; empty ones stay null
		B b(wala::with_capacity{2}, 2, tracked("x"));
		const tracked* p = b.data();
		V v = std::move(b).into_vec();
		REQUIRE(v.data() == p);
		REQUIRE(show(v) == "x,x,");
		S s(wala::with_capacity{2}, 2, tracked("y"));
		p = s.data();
		V v2 = std::move(s).into_vec();
		REQUIRE(v2.data() == p);
		REQUIRE(show(v2) == "y,y,");
		V e;
		B eb = std::move(e).into_bounded();
		REQUIRE(eb.data() == nullptr);
		REQUIRE(eb.capacity() == 0);
		G eg = std::move(eb).into_growable();
		eg.emplace_back("z");
		REQUIRE(show(eg) == "z,");
	}
	REQUIRE(tracked::live == 0);

	wala::bounded_vec<std::unique_ptr<int>> p(wala::with_capacity{2});
	p.emplace_back(std::make_unique<int>(1));
	wala::growable_vec<std::unique_ptr<int>> g = std::move(p).into_growable();
	g.emplace_back(std::make_unique<int>(2));
	wala::vec<std::unique_ptr<int>> v = std::move(g).into_vec();
	REQUIRE(v.size() == 2);
	REQUIRE(*v[1] == 2);

	require_debug_abort([] { wala::bounded_vec<int> b(wala::with_capacity{3}, 2); (void)std::move(b).into_vec(); });
	require_debug_abort([] { wala::bounded_stack<int> b(wala::with_capacity{3}, 2); (void)std::move(b).into_vec(); });
}
