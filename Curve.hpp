#pragma once

#include <ostream>
#include <vector>
#include <initializer_list>
#include <utility>

struct Point {
	double x, y, rate; // in kbits per millisecond = mbits per second

	Point(double, double, double) noexcept;

	Point operator-() const noexcept;
};

class Curve;
extern Curve max_curve(const Curve&, const Curve&);
extern Curve min_curve(const Curve&, const Curve&);
extern double max_h_distance(const Curve&, const Curve&) noexcept;
extern double max_v_distance(const Curve&, const Curve&) noexcept;

class Curve: protected std::vector<Point> {
private:
	typedef std::vector<Point> base_type;

public:
	typedef base_type::iterator               iterator;
	typedef base_type::const_iterator         const_iterator;
	typedef base_type::reverse_iterator       reverse_iterator;
	typedef base_type::const_reverse_iterator const_reverse_iterator;

	typedef base_type::reference              reference;
	typedef base_type::const_reference        const_reference;

	typedef base_type::pointer                pointer;
	typedef base_type::const_pointer          const_pointer;

	typedef base_type::value_type             value_type;
	typedef base_type::size_type              size_type;
	typedef base_type::difference_type        difference_type;

	Curve() noexcept = default;
	explicit Curve(Point&&);
	Curve(Point&&, Point&&);
	Curve(std::initializer_list<Point>);

	Curve(const Curve&)                = default;
	Curve(Curve&&) noexcept            = default;
	Curve& operator=(const Curve&)     = default;
	Curve& operator=(Curve&&) noexcept = default;

	Curve operator-() const;

	Curve operator+(const Point&) const;
	Curve operator-(const Point&) const;
	Curve& operator+=(const Point&);
	Curve& operator-=(const Point&);

	Curve operator+(const Curve&) const;
	Curve operator-(const Curve&) const;

	Curve& operator+=(const Curve&);
	Curve& operator-=(const Curve&);

	using base_type::operator[];
	using base_type::front;
	using base_type::back;
	using base_type::begin;
	using base_type::end;
	using base_type::cbegin;
	using base_type::cend;
	using base_type::rbegin;
	using base_type::rend;
	using base_type::crbegin;
	using base_type::crend;
	using base_type::emplace_back;

	using base_type::size;

	double at(double x) const noexcept;

	double x(double y) const noexcept;

private:
	Curve(base_type&&) noexcept;

	Curve& merge(const Point&);
	Curve& merge(const Curve&);
	Curve& merge_neg(const Curve&);
	unsigned int merge_helper(const Point&);

friend Curve max_curve(const Curve&, const Curve&);
friend Curve min_curve(const Curve&, const Curve&);
// friend double max_h_distance(const Curve&, const Curve&);
// friend double max_v_distance(const Curve&, const Curve&);
};
extern std::ostream& operator<<(std::ostream& os, const Curve&);


inline
Point::Point(double x__, double y__, double rate__) noexcept
	: x(x__)
	, y(y__)
	, rate(rate__)
{
}

inline
Curve::Curve(Point&& p)
	: base_type{std::move(p)}
{
}

inline
Curve::Curve(Point&& p1, Point&& p2)
	: base_type{std::move(p1), std::move(p2)}
{
}

inline
Curve::Curve(std::initializer_list<Point> l)
	: base_type(l)
{
}

inline
Curve::Curve(base_type&& v) noexcept
	: base_type(std::move(v))
{
}

inline Point
Point::operator-() const noexcept {
	return Point(x, -y, -rate);
}
