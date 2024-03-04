#include "Curve.hpp"
#include <utility>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cassert>

static constexpr double PSEUDO_ZERO = 1e-7;
typedef std::vector<Point> Points;

static bool max_curve1(Curve::const_iterator& i, Curve::const_iterator __li,
		       Curve::const_iterator& j, Curve::const_iterator __lj,
		       Points& r);
static bool max_curve2(Curve::const_iterator& i, Curve::const_iterator __li,
		       Curve::const_iterator& j, Curve::const_iterator __lj,
		       Points& r);
static bool max_curve_tail(Curve::const_iterator& i, Points& r);
static bool max_curve3(const Point&, const Point&, Points& r);
static bool max_curve4(const Point&, const Point&, Points& r);

static bool min_curve1(Curve::const_iterator& i, Curve::const_iterator __li,
		       Curve::const_iterator& j, Curve::const_iterator __lj,
		       Points& r);
static bool min_curve2(Curve::const_iterator& i, Curve::const_iterator __li,
		       Curve::const_iterator& j, Curve::const_iterator __lj,
		       Points& r);
static bool min_curve_tail(Curve::const_iterator& i, Points& r);
static bool min_curve3(const Point&, const Point&, Points& r);
static bool min_curve4(const Point&, const Point&, Points& r);
static Point intersect(const Point& p1, const Point& p2) noexcept;

Curve&
Curve::operator+=(const Point& p) {
	if (!empty())
		return merge(p);
	push_back(p);
	return *this;
}

Curve&
Curve::operator+=(const Curve& c) {
	if (!c.empty()) {
		if (c.size() > 1) {
			if (!empty())
				return merge(c);
			return *this = c;
		}
		return *this += c[0];
	}
	return *this;
}

Curve
Curve::operator+(const Point& p) const {
	Curve r(*this);
	return r += p;
}

Curve
Curve::operator+(const Curve& c) const {
	Curve r(*this);
	return r += c;
}

Curve&
Curve::operator-=(const Point& p) {
	if (!empty())
		return merge(-p);
	push_back(-p);
	return *this;
}

Curve&
Curve::operator-=(const Curve& c) {
	if (!c.empty()) {
		if (c.size() > 1) {
			if (!empty())
				return merge_neg(c);
			reserve(c.size());
			for (const Point& p: c)
				push_back(-p);
			return *this;
		}
		return *this -= c[0];
	}
	return *this;
}

Curve
Curve::operator-(const Point& p) const {
	Curve r(*this);
	return r -= p;
}

Curve
Curve::operator-(const Curve& c) const {
	Curve r(*this);
	return r -= c;
}

Curve
Curve::operator-() const {
	Curve r(*this);
	for (Point& p: r) {
		p = -p;
	}
	return r;
}

double
Curve::y(double x) const noexcept {
	const_iterator i = std::lower_bound(begin(), end(), x,
					    [](const Point& p, double x) {
						    return p.x < x;
					    });
	if (end() == i || i->x != x) {
		if (begin() != i) {
			--i;
			return i->y + (x - i->x) * i->rate;
		}
		return 0.0;
	}
	return i->y;
}

double
Curve::x(double y) const noexcept {
	unsigned int j = 0, L = size();
	while (j < L) {
		if ((*this)[j].rate != 0.0) {
			double u = (*this)[j].x + (y - (*this)[j].y) / (*this)[j].rate;
			if (u >= (*this)[j++].x)
				if (j >= L || u < (*this)[j].x)
					return u;
		} else {
			if ((*this)[j].y == y)
				return (*this)[j].x;
			++j;
		}
	}
	return 0.0;
}

Curve
max_curve(const Curve& c1, const Curve& c2) {
	if (c1.empty())
		return c2;
	if (c2.empty())
		return c1;

	Points r;

	Curve::const_iterator i = c1.begin(), __li = c1.end(),
		j = c2.begin(), __lj = c2.end();

	bool added_from_c1 = max_curve1(i, __li, j, __lj, r);

	while (i != __li && j != __lj) {
		if (i->x < j->x) {
			if (added_from_c1)
				// c1 c1
				added_from_c1 = max_curve3(*i, *(j - 1), r);
			else
				// c2 c1
				added_from_c1 = max_curve4(*i, *(i - 1), r);
			++i;
		} else if (j->x < i->x) {
			if (added_from_c1)
				// c1 c2
				added_from_c1 = !max_curve4(*j, *(j - 1), r);
			else
				// c2 c2
				added_from_c1 = !max_curve3(*j, *(i - 1), r);
			++j;
		} else {
			// we assume that the curves are continuous
			if (i->y < j->y || (i->y == j->y && i->rate < j->rate)) {
				r.push_back(*j);
				added_from_c1 = false;
			} else {
				r.push_back(*i);
				added_from_c1 = true;
			}
			++i;
			++j;
		}
	}

	for (; __li != i; ++i)
		if (added_from_c1)
			// c1 c1
			added_from_c1 = max_curve3(*i, *(j - 1), r);
		else
			// c2 c1
			added_from_c1 = max_curve4(*i, *(i - 1), r);

	for (; __lj != j; ++j)
		if (added_from_c1)
			// c1 c2
			added_from_c1 = !max_curve4(*j, *(j - 1), r);
		else
			// c2 c2
			added_from_c1 = !max_curve3(*j, *(i - 1), r);

	const Point& lastI = *(i - 1);
	const Point& lastJ = *(j - 1);
	if (added_from_c1) {
		if (lastJ.rate > lastI.rate)
			r.push_back(intersect(lastJ, lastI));
	} else {
		if (lastI.rate > lastJ.rate)
			r.push_back(intersect(lastI, lastJ));
	}
	return Curve(std::move(r));
}

// pre-conditions:
// i == c1.begin()	i is a valid iterator and points to the first point of c1
// i != c1.end()
// j == c2.begin()	j is a valid iterator and points to the first point of c2
// j != c2.end()
// r.empty()		the result curve has no points

// post-conditions:
// i != c1.begin() && j == c2.begin() + 1 || i == c1.begin() + 1 && j != c2.begin()
static bool
max_curve1(Curve::const_iterator& i, Curve::const_iterator __li,
	   Curve::const_iterator& j, Curve::const_iterator __lj, Points& r) {
	if (i->x < j->x) {
		r.push_back(*i);
		++i;
		if (__li == i)
			return max_curve_tail(j, r);
		return max_curve2(i, __li, j, __lj, r);
	}
	if (i->x > j->x) {
		r.push_back(*j);
		++j;
		if (__lj == j)
			return !max_curve_tail(i, r);
		return !max_curve2(j, __lj, i, __li, r);
	}
	// i->x == j->x
	if (i->y < j->y || (i->y == j->y && i->rate < j->rate)) {
		r.push_back(*j);
		++i;
		++j;
		return false;
	}
	// i->y > j->y || i->y == j->y && i->rate >= j->rate
	r.push_back(*i);
	++i;
	++j;
	return true;
}

// pre-conditions:
// r.size() == 1	the output curve contains exactly one point
// i == c1.begin() + 1	i is a valid iterator and points to the second point of c1
// i != __li
// j == c2.begin()	j is a valid iterator and points to the first point of c2
// j != __lj

// post-conditions:
// j == c2.begin() + 1
static bool
max_curve2(Curve::const_iterator& i, Curve::const_iterator __li,
	   Curve::const_iterator& j, Curve::const_iterator __lj, Points& r) {
	while (i->x < j->x) {
		r.push_back(*i);
		++i;
		if (__li == i)
			return max_curve_tail(j, r);
	}
	if (i->x > j->x)
		return max_curve_tail(j, r);
	// i->x == j->x
	if (i->y < j->y || (i->y == j->y && i->rate < j->rate)) {
		r.push_back(*j);
		++i;
		++j;
		return false;
	}
	// i->x == j->x
	// i->y > j->y || i->y == j->y && i->rate >= j->rate
	r.push_back(*i);
	++i;
	++j;
	return true;
}

// pre-conditions:
// i == c1.begin()	i is a valid iterator and points to the first point of c1
// i != c1.end()
// !r.empty()		r contains at least one point

// post-conditions
// i == c1.begin() + 1
static bool
max_curve_tail(Curve::const_iterator& i, Points& r) {
	const Point& prev = r.back();
	double y = prev.y + prev.rate * (i->x - prev.x);
	if (y < i->y || (y == i->y && prev.rate < i->rate)) {
		r.push_back(*i);
		++i;
		return false;
	}
	// y > i->y || y == i->y && prev.rate >= i->rate
	++i;
	return true;
}

// the last point of r comes from curve1
// I is on curve 1
// J is on curve 2
static bool
max_curve3(const Point& I, const Point& J, Points& r) {
	// curve1 curve1
	double y = J.y + J.rate * (I.x - J.x);
	if (y > I.y) {
		// c2 overtakes
		r.push_back(intersect(J, r.back()));
		return false;
	}
	if (y == I.y && J.rate > I.rate) {
		// c2 equalizes. Its rate is higher than
		// the new rate of c1 at I.x.
		r.emplace_back(I.x, I.y, J.rate);
		return false;
	}
	r.push_back(I);
	return true;
}

// the last point of r comes from curve2
// I is on curve1
// Iprev is on curve1
static bool
max_curve4(const Point& I, const Point& Iprev, Points& r) {
	// c2 c1
	const Point& J = r.back();
	double y = J.y + J.rate * (I.x - J.x);
	if (I.y > y) {
		// c1 overtakes c2
		r.push_back(intersect(Iprev, J));
		r.push_back(I);
		return true;
	}
	if (I.y == y && I.rate > J.rate) {
		// c1 equalizes. Its rate is higher than
		// the rate of c2 at I.x.
		r.push_back(I);
		return true;
	}
	// else c1 does not overtake c2 at this point
	return false;
}

Curve
min_curve(const Curve& c1, const Curve& c2) {
	if (c1.empty())
		return c2;
	if (c2.empty())
		return c1;

	Points r;

	Curve::const_iterator i = c1.begin(), __li = c1.end(),
		j = c2.begin(), __lj = c2.end();

	bool added_from_c1 = min_curve1(i, __li, j, __lj, r);

	while (i != __li && j != __lj) {
		if (i->x < j->x) {
			if (added_from_c1)
				// c1 c1
				added_from_c1 = min_curve3(*i, *(j - 1), r);
			else
				// c2 c1
				added_from_c1 = min_curve4(*i, *(i - 1), r);
			++i;
		} else if (j->x < i->x) {
			if (added_from_c1)
				// c1 c2
				added_from_c1 = !min_curve4(*j, *(j - 1), r);
			else
				// c2 c2
				added_from_c1 = !min_curve3(*j, *(i - 1), r);
			++j;
		} else {
			// we assume that the curves are continuous
			if (j->y < i->y || (i->y == j->y && j->rate < i->rate)) {
				r.push_back(*j);
				added_from_c1 = false;
			} else {
				r.push_back(*i);
				added_from_c1 = true;
			}
			++i;
			++j;
		}
	}

	for (; __li != i; ++i)
		if (added_from_c1)
			// c1 c1
			added_from_c1 = min_curve3(*i, *(j - 1), r);
		else
			// c2 c1
			added_from_c1 = min_curve4(*i, *(i - 1), r);

	for (; __lj != j; ++j)
		if (added_from_c1)
			// c1 c2
			added_from_c1 = !min_curve4(*j, *(j - 1), r);
		else
			// c2 c2
			added_from_c1 = !min_curve3(*j, *(i - 1), r);

	const Point& lastI = *(i - 1);
	const Point& lastJ = *(j - 1);
	if (added_from_c1) {
		if (lastJ.rate < lastI.rate)
			r.push_back(intersect(lastJ, lastI));
	} else {
		if (lastI.rate < lastJ.rate)
			r.push_back(intersect(lastI, lastJ));
	}

	return Curve(std::move(r));
}

// pre-conditions:
// i == c1.begin()	i is a valid iterator and points to the first point of c1
// i != c1.end()
// j == c2.begin()	j is a valid iterator and points to the first point of c2
// j != c2.end()
// r.empty()		the result curve has no points

// post-conditions:
// i != c1.begin() && j == c2.begin() + 1 || i == c1.begin() + 1 && j != c2.begin()
static bool
min_curve1(Curve::const_iterator& i, Curve::const_iterator __li,
	   Curve::const_iterator& j, Curve::const_iterator __lj, Points& r) {
	if (i->x < j->x) {
		r.push_back(*i);
		++i;
		if (__li == i)
			return min_curve_tail(j, r);
		return min_curve2(i, __li, j, __lj, r);
	}
	if (i->x > j->x) {
		r.push_back(*j);
		++j;
		if (__lj == j)
			return !min_curve_tail(i, r);
		return !min_curve2(j, __lj, i, __li, r);
	}
	// i->x == j->x
	if (j->y < i->y || (j->y == i->y && j->rate < i->rate)) {
		r.push_back(*j);
		++i;
		++j;
		return false;
	}
	// i->y < j->y || i->y == j->y && i->rate <= j->rate
	r.push_back(*i);
	++i;
	++j;
	return true;
}

// pre-conditions:
// r.size() == 1	the output curve contains exactly one point
// i == c1.begin() + 1	i is a valid iterator and points to the second point of c1
// i != __li
// j == c2.begin()	j is a valid iterator and points to the first point of c2
// j != __lj

// post-conditions:
// j == c2.begin() + 1
static bool
min_curve2(Curve::const_iterator& i, Curve::const_iterator __li,
	   Curve::const_iterator& j, Curve::const_iterator __lj, Points& r) {
	while (i->x < j->x) {
		r.push_back(*i);
		++i;
		if (__li == i)
			return min_curve_tail(j, r);
	}
	if (i->x > j->x)
		return min_curve_tail(j, r);
	// i->x == j->x
	if (j->y < i->y || (j->y == i->y && j->rate < i->rate)) {
		r.push_back(*j);
		++i;
		++j;
		return false;
	}
	// i->x == j->x
	// i->y < j->y || i->y == j->y && i->rate <= j->rate
	r.push_back(*i);
	++i;
	++j;
	return true;
}

// pre-conditions:
// i == c1.begin()	i is a valid iterator and points to the first point of c1
// i != c1.end()
// !r.empty()		r contains at least one point

// post-conditions
// i == c1.begin() + 1
static bool
min_curve_tail(Curve::const_iterator& i, Points& r) {
	const Point& prev = r.back();
	double y = prev.y + prev.rate * (i->x - prev.x);
	if (y > i->y || (y == i->y && prev.rate > i->rate)) {
		r.push_back(*i);
		++i;
		return false;
	}
	// y < i->y || y == i->y && prev.rate <= i->rate
	++i;
	return true;
}

// the last point of r comes from curve1
// I is on curve 1
// J is on curve 2
static bool
min_curve3(const Point& I, const Point& J, Points& r) {
	// curve1 curve1
	double y = J.y + J.rate * (I.x - J.x);
	if (y < I.y) {
		// c2 overtakes
		r.push_back(intersect(J, r.back()));
		return false;
	}
	if (y == I.y && J.rate < I.rate) {
		// c2 equalizes. Its rate is lower than
		// the new rate of c1 at I.x.
		r.emplace_back(I.x, I.y, J.rate);
		return false;
	}
	r.push_back(I);
	return true;
}

// the last point of r comes from curve2
// I is on curve1
// Iprev is on curve1
static bool
min_curve4(const Point& I, const Point& Iprev, Points& r) {
	// c2 c1
	const Point& J = r.back();
	double y = J.y + J.rate * (I.x - J.x);
	if (I.y < y) {
		// c1 overtakes c2
		r.push_back(intersect(Iprev, J));
		r.push_back(I);
		return true;
	}
	if (I.y == y && I.rate < J.rate) {
		// c1 equalizes. Its rate is lower than
		// the rate of c2 at I.x.
		r.push_back(I);
		return true;
	}
	// else c1 does not overtake c2 at this point
	return false;
}

double
max_h_distance(const Curve& f_, const Curve& g_) noexcept {
	assert(f_.size() > 0);
	assert(g_.size() > 0);

	const Point *f = &f_.back(), *fb = &f_.front();
	assert(f->rate >= 0.0);
	const Point *g = &g_.back(), *gb = &g_.front();
	assert(g->rate >= 0.0);

	if (f->rate > g->rate || (0.0 == g->rate && f->y > g->y))
		return std::numeric_limits<double>::infinity();

	const Point *fnxt = nullptr, *gnxt = nullptr;
	double r = 0.0;
	double dx = g->x - f->x;
	double dy = g->y - f->y;
	// f->rate <= g->rate && g->rate > 0.0 ||
	// f->rate == g->rate == 0.0 && f->y <= g->y
	do {
		double dx = g->x - f->x;
		double dy = g->y - f->y;
		if (dy >= 0.0) { // g.y >= f.y
			if (f->rate <= g->rate && f->rate > 0.0) {
				double q = dx - dy / f->rate;
				if (q > r)
					r = q;
			}
			if (gb == g) {
				double df = gb->x - fb->x;
				return df > r ? df : r;
			}
			gnxt = g--;
			assert(g->rate > gnxt->rate);
			assert(fabs(g->rate * (gnxt->x - g->x) + g->y - gnxt->y) < PSEUDO_ZERO);
		} else { // f.y > g.y
			if (f->rate <= g->rate) {
				double q = dx - dy / g->rate;
				if (q > r)
					r = q;
			}
			if (fb == f)
				return r;
			fnxt = f--;
			assert(f->rate > fnxt->rate);
			assert(fabs(f->rate * (fnxt->x - f->x) + f->y - fnxt->y) < PSEUDO_ZERO);
		}
	} while (true);
}

double
max_v_distance(const Curve& f_, const Curve& g_) noexcept {
	assert(f_.size() > 0);
	assert(g_.size() > 0);

	const Point *f = &f_.back(), *fb = &f_.front();
	assert(f->rate >= 0.0);
	const Point *g = &g_.back(), *gb = &g_.front();
	assert(g->rate >= 0.0);

	if (f->rate > g->rate)
		return std::numeric_limits<double>::infinity();

	const Point *fnxt = nullptr, *gnxt = nullptr;
	double r = 0.0;
	// f->rate <= g->rate

	do {
		double dx = g->x - f->x;
		double dy = g->y - f->y;
		if (dx >= 0.0) { // g.x >= f.x
			if (f->rate <= g->rate) {
				double q = -dy + dx * f->rate;
				if (q > r)
					r = q;
			}
			if (gb == g) {
				if (0.0 != dx || fb != f) {
					double q = f->y + dx * f->rate;
					if (q > r)
						return q;
				}
				return r;
			}
			gnxt = g--;
			assert(g->rate > gnxt->rate);
			assert(fabs(g->rate * (gnxt->x - g->x) + g->y - gnxt->y) < PSEUDO_ZERO);
		} else { // f.x > g.x
			if (f->rate <= g->rate) {
				double q = -dy + dx * g->rate;
				if (q > r)
					r = q;
			}
			if (fb == f)
				return r;
			fnxt = f--;
			assert(f->rate > fnxt->rate);
			assert(fabs(f->rate * (fnxt->x - f->x) + f->y - fnxt->y) < PSEUDO_ZERO);
		}
	} while (true);
	return r;
}

Curve
max_conv(const Curve& c1, const Curve& c2) {
	// (f x g)(t) \sup_{0 \le \lambda \le t} { f(t - \lambda) + g(\lambda) }
	return Curve{};
}

Curve
max_deconv(const Curve& c1, const Curve& c2) {
	// (f x g)(t) \inf_{\lambda \ge 0} { f(t + \lambda) - g(\lambda) }
	return Curve{};
}

Curve
min_conv(const Curve& c1, const Curve& c2) {
	// (f x g)(t) \inf_{0 \le \lambda \le t} { f(t - \lambda) + g(\lambda) }
	return Curve{};
}

Curve
min_deconv(const Curve& c1, const Curve& c2) {
	// (f x g)(t) \sup_{\lambda \ge 0} { f(t + \lambda) - g(\lambda) }
	return Curve{};
}

static Point
intersect(const Point& p1, const Point& p2) noexcept {
	return Point((p2.y - p1.y + p1.rate * p1.x - p2.rate * p2.x) / (p1.rate - p2.rate),
		     (p1.rate * p2.y - p2.rate * p1.y + p1.rate * p2.rate * (p1.x - p2.x)) / (p1.rate - p2.rate),
		     p1.rate);
}

Curve&
Curve::merge(const Point& p) {
	unsigned int j = merge_helper(p);
	unsigned int L = size();
	for (; j < L; ++j) {
		(*this)[j].y    += p.y + p.rate * ((*this)[j].x - p.x);
		(*this)[j].rate += p.rate;
	}
	return *this;
}

unsigned int
Curve::merge_helper(const Point& p) {
	iterator i = std::lower_bound(begin(), end(), p.x,
				      [](const Point& p, double x) {
					      return p.x < x;
				      });
	unsigned int r = i - begin() + 1;
	if (end() == i || i->x != p.x) {
		if (begin() != i) {
			--i;
			Point n(p.x,
				i->y + p.y + i->rate * (p.x - i->x),
				i->rate + p.rate);
			insert(++i, std::move(n));
		} else
			insert(i, p);
	} else {
		i->y    += p.y;
		i->rate += p.rate;
	}
	return r;
}

Curve&
Curve::merge(const Curve& c) {
	// neither this nor c are empty
	// merge_helper merges the first point of c into this
	// j cannot be 0 after this call
	unsigned int j = merge_helper(c[0]);
	// j > 0

	// (A): c[i-1].x <= this[j-1].x
	// (B): c[i-1].x = this[j-1].x
	// (obviously (B) implies (A))
	unsigned int i = 1, L = size(), I = c.size();
	for (; j < L && i < I; ++j) {
		if ((*this)[j].x < c[i].x) {
			// this[j-1].x < this[j].x
			// (A)
			// =======================
			// => c[i-1].x <= this[j-1].x < this[j].x => this[j].x - c[i-1].x > 0
			(*this)[j].y    += c[i-1].y + c[i-1].rate * ((*this)[j].x - c[i-1].x);
			(*this)[j].rate += c[i-1].rate;
			// j is incremented in the for loop
			// (A) stays true
		} else if ((*this)[j].x > c[i].x) {
			// c[i-1].x < c[i].x
			// (B)
			// ==================
			// this[j-1].x < c[i].x => c[i].x - this[j-1].x > 0
			Point p(c[i].x,
				(*this)[j-1].y + (*this)[j-1].rate * (c[i].x - (*this)[j-1].x),
				(*this)[j-1].rate + c[i].rate - c[i-1].rate);
			insert(begin() + j, std::move(p));
			++L;
			++i;
			// j is incremented in the for loop
			// (A) stays true
			// (B) remains true
		} else {
			(*this)[j].y    += c[i].y;
			(*this)[j].rate += c[i].rate;
			++i;
			// j is incremented in the for loop
			// (A) stays true
			// (B) stays/becomes true
		}
	}
	for (; j < L; ++j) {
		(*this)[j].y    += c[i-1].y + c[i-1].rate * ((*this)[j].x - c[i-1].x);
		(*this)[j].rate += c[i-1].rate;
	}
	// j == L
	// i > 0
	for (; i < I; ++i, ++j) {
		Point p(c[i].x,
			(*this)[j-1].y + (*this)[j-1].rate * (c[i].x - (*this)[j-1].x),
			(*this)[j-1].rate + c[i].rate - c[i-1].rate);
		emplace_back(std::move(p));
	}
	return *this;
}

Curve&
Curve::merge_neg(const Curve& c) {
	// neither this nor c are empty
	// merge_helper merges the first point of c into this
	// j cannot be 0 after this call
	unsigned int j = merge_helper(-c[0]);
	// j > 0

	// (A): c[i-1].x <= this[j-1].x
	// (B): c[i-1].x = this[j-1].x
	// (obviously (B) implies (A))
	unsigned int i = 1, L = size(), I = c.size();
	for (; j < L && i < I; ++j) {
		if ((*this)[j].x < c[i].x) {
			// this[j-1].x < this[j].x
			// (A)
			// =======================
			// => c[i-1].x <= this[j-1].x < this[j].x => this[j].x - c[i-1].x > 0
			(*this)[j].y    -= c[i-1].y + c[i-1].rate * ((*this)[j].x - c[i-1].x);
			(*this)[j].rate -= c[i-1].rate;
			// j is incremented in the for loop
			// (A) stays true
		} else if ((*this)[j].x > c[i].x) {
			// c[i-1].x < c[i].x
			// (B)
			// ==================
			// this[j-1].x < c[i].x => c[i].x - this[j-1].x > 0
			Point p(c[i].x,
				(*this)[j-1].y + (*this)[j-1].rate * (c[i].x - (*this)[j-1].x),
				(*this)[j-1].rate - c[i].rate + c[i-1].rate);
			insert(begin() + j, std::move(p));
			++L;
			++i;
			// j is incremented in the for loop
			// (A) stays true
			// (B) remains true
		} else {
			(*this)[j].y    -= c[i].y;
			(*this)[j].rate -= c[i].rate;
			++i;
			// j is incremented in the for loop
			// (A) stays true
			// (B) stays/becomes true
		}
	}
	for (; j < L; ++j) {
		(*this)[j].y    -= c[i-1].y + c[i-1].rate * ((*this)[j].x - c[i-1].x);
		(*this)[j].rate -= c[i-1].rate;
	}
	// j == L
	// i > 0
	for (; i < I; ++i, ++j) {
		Point p(c[i].x,
			(*this)[j-1].y + (*this)[j-1].rate * (c[i].x - (*this)[j-1].x),
			(*this)[j-1].rate - c[i].rate + c[i-1].rate);
		emplace_back(std::move(p));
	}
	return *this;
}

// for gnuplot
std::ostream&
operator<<(std::ostream& os, const Curve& c) {
	Curve::const_iterator j = c.begin(), __li = c.end();
	if (__li != j) {
		Curve::const_iterator i = j++;
		os << "x < " << i->x << " ? 0 : ";
		if (__li != j) {
			os << "x < " << j->x << " ? " << i->y << "+(x-" << i->x << ")*" << i->rate << " : ";
			for (++i, ++j; __li != j; ++i, ++j)
				os << "x < " << j->x << " ? " << i->y << "+(x-" << i->x << ")*" << i->rate << " : ";
		}
		os << i->y << "+(x-" << i->x << ")*" << i->rate;
	}
	return os;
}
