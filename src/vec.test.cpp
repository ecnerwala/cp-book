#include <catch2/catch_test_macros.hpp>

#include "vec.hpp"

#include <string>

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

TEST_CASE("vec", "[wala::vec]") {
	// TODO
}
