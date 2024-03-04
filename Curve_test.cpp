#include "Curve.hpp"
#include <iostream>
#include <fstream>
#include <random>
#include <ctime>
#include <cassert>

template<typename... C>
static Point recalculate(C const&... curves);
static Curve random_curve(unsigned int segments, double len);

static std::mt19937_64 g;

int
main(int argc, char *argv[]) {
	g.seed(time(nullptr));
	// Curve c1{
	// 	Point{ 0.0,   1.0, 28.0},
	// 	Point{ 1.0,  29.0, 26.0},
	// 	Point{ 1.5,  42.0, 24.0},
	// 	Point{ 3.0,  78.0, 22.0},
	// 	Point{ 8.0, 188.0, 20.0},
	// 	Point{11.5, 258.0, 18.0},
	// 	Point{15.0, 321.0, 16.0},
	// 	Point{19.0, 385.0, 14.0},
	// 	Point{20.0, 399.0, 12.0},
	// 	Point{25.0, 459.0, 10.0},
	// 	Point{27.0, 479.0,  8.0},
	// 	Point{29.0, 495.0,  6.0},
	// 	Point{31.0, 507.0,  4.0},
	// 	Point{34.0, 519.0,  2.0}
	// };
	// Curve c2{
	// 	Point{ 0.0,   3.0, 14.0},
	// 	Point{ 4.0,  59.0, 13.0},
	// 	Point{ 5.0,  72.0, 12.0},
	// 	Point{ 7.0,  96.0, 11.0},
	// 	Point{10.0, 129.0, 10.0},
	// 	Point{13.5, 164.0,  9.0},
	// 	Point{18.0, 204.5,  8.0},
	// 	Point{22.0, 236.5,  7.0},
	// 	Point{24.0, 250.5,  6.0},
	// 	Point{25.0, 256.5,  5.0},
	// 	Point{28.0, 271.5,  4.0}
	// };
	// Curve c3{
	// 	Point{ 8.0,   1.0, 28.0},
	// 	Point{ 9.0,  29.0, 26.0},
	// 	Point{ 9.5,  42.0, 24.0},
	// 	Point{11.0,  78.0, 22.0},
	// 	Point{16.0, 188.0, 20.0},
	// 	Point{19.5, 258.0, 18.0},
	// 	Point{23.0, 321.0, 16.0},
	// 	Point{27.0, 385.0, 14.0},
	// 	Point{28.0, 399.0, 12.0},
	// 	Point{33.0, 459.0, 10.0},
	// 	Point{35.0, 479.0,  8.0},
	// 	Point{37.0, 495.0,  6.0},
	// 	Point{39.0, 507.0,  4.0},
	// 	Point{42.0, 519.0,  2.0}
	// };

	// std::cout << c1 << std::endl;
	// std::cout << c2 << std::endl;
	// std::cout << c3 << std::endl;

	// std::cout << c1 + c2 << std::endl;
	// std::cout << c2 + c1 << std::endl;
	// std::cout << c3 + c2 << std::endl;
	// std::cout << c2 + c3 << std::endl;

	// Curve diff = c1 - c2;
	// std::cout << diff << std::endl;

	// std::ofstream f("plot");
	// for (double x = 0.0; x < 34.0; x += 0.05) {
	// 	f << x << '\t' << diff.y(x) << std::endl;
	// }
	// f.close();

	// Point p1 = recalculate(c1, c2, c3);
	// const Curve& r = c1 + c2 + c3;
	// const Point& p2 = r.back();

	// std::cout << p1.x << '\t' << p1.y << '\t' << p1.rate << std::endl;
	// std::cout << p2.x << '\t' << p2.y << '\t' << p2.rate << std::endl;

	// Curve c1{
	// 	Point{0, 2, 2.5},
	// 	Point{2, 7, 1},
	// 	Point{5, 10, 5.0/6},
	// 	Point{11, 15, 4.0/9},
	// 	Point{20, 19, 1.0/3}
	// };

	// Curve c2{
	// 	Point{0, 5, 3.5/3},
	// 	Point{3, 8.5, 2.7/3},
	// 	Point{6, 11.2, 4.3/7},
	// 	Point{13, 15.5, 0.575},
	// 	Point{17, 17.8, 1.6/3},
	// 	Point{23, 21, 1.0/5}
	// };

	// std::cout << c1 << std::endl;
	// std::cout << c2 << std::endl;

	// Curve c3 = max_curve(c1, c2);
	// std::cout << c3 << std::endl;

	// Curve c4 = min_curve(c1, c2);
	// std::cout << c4 << std::endl;

	// Curve c1{
	// 	Point{                     0.0,                     50.0, 5.67128181961770953035}, // 80
	// 	Point{ 17.36481776669303488700, 148.48077530122080593600, 2.74747741945462227867}, // 70
	// 	Point{ 51.56683209925990819200, 242.45003737981164434000, 1.73205080756887729352}, // 60
	// 	Point{101.56683209925990819200, 329.05257775825550901600, 1.19175359259420995867}, // 50
	// 	Point{165.84559306791384082500, 405.65702207015331253500, 0.83909963117728001174}, // 40
	// 	Point{242.45003737981164434500, 469.93578303880724516600, 0.57735026918962576451}, // 30
	// 	Point{329.05257775825550902100, 519.93578303880724516600, 0.36397023426620236133}, // 20
	// 	Point{423.02183983684634742600, 554.13779737137411846900, 0.17632698070846497346}  // 10
	// };
	// Curve c2{
	// 	Point{ -2.18304714542560309999,                      0.0, 11.43005230276134306558}, // 85
	// 	Point{  6.53252712934021425701,  99.61946980917455322900,  3.73205080756887729332}, // 75
	// 	Point{ 32.41443163959229049301, 196.21205243808138190300,  2.14450692050955861631}, // 65
	// 	Point{ 74.67625781366223411201, 286.84283114174637822600,  1.42814800674211450213}, // 55
	// 	Point{132.03390144876684372301, 368.75803557064555719300,  1.0},                    // 45
	// 	Point{202.74457956742159616301, 439.46871368930030963300,  0.70020753820970977943}, // 35
	// 	Point{284.65978399632077513201, 496.82635732440491924200,  0.46630765815499859281}, // 25
	// 	Point{375.29056269998577145601, 539.08818349847486285900,  0.26794919243112270646}, // 15
	// };

	// std::cout << c1 << std::endl;
	// std::cout << c2 << std::endl;
	// std::cout << max_h_distance(c1, c2) << std::endl;

	Curve c1 = random_curve(8, 100);
	Curve c2 = random_curve(8, 100);
	const Curve *cp1 = &c1, *cp2 = &c2;;
	if (c1.back().rate > c2.back().rate) {
		cp1 = &c2;
		cp2 = &c1;
	}
	/*
	  x = (f.x * f.rate - g.x * g.rate - f.y + g.y) / (f.rate - g.rate)
	  y = (g.y * f.rate - f.y * g.rate + (f.x - g.x) * f.rate * g.rate) / (f.rate - g.rate)
	 */
	const Point& p1 = cp1->back();
	const Point& p2 = cp2->back();
	double hi = 0.0;
	if (p1.rate != p2.rate) {
		double x = (p1.x * p1.rate - p2.x * p2.rate - p1.y + p2.y) / (p1.rate - p2.rate);
		double y = (p2.y * p1.rate - p1.y * p2.rate + (p1.x - p2.x) * p1.rate * p2.rate) /
			(p1.rate - p2.rate);
		if (x >= p1.x && x >= p2.x)
			hi = y;
		else
			hi = p1.y > p2.y ? p1.y : p2.y;
	} else
		hi = p1.y > p2.y ? p1.y : p2.y;
	std::cout << max_h_distance(*cp1, *cp2) << std::endl;
	double r = 0.0;
	for (double y = 0.0; y <= hi; y += 1e-4) {
		double diff = cp2->x(y) - cp1->x(y);
		if (diff > r)
			r = diff;
	}
	std::cout << r << std::endl;

	r = -std::numeric_limits<double>::max();
	std::cout << max_v_distance(*cp1, *cp2) << std::endl;
	hi = p1.x > p2.x ? p1.x + 10.0 : p2.x + 10.0;
	for (double x = 0.0; x <= hi; x += 1e-4) {
		double y2 = cp2->y(x);
		double y1 = cp1->y(x);
		std::cerr << x << '\t' << y1 << '\t' << y2 << std::endl;
		double diff = y1 - y2;
		if (diff > r)
			r = diff;
	}
	std::cout << r << std::endl;

	std::cout << c1 << std::endl;
	std::cout << c2 << std::endl;

	return 0;
}

static Curve
random_curve(unsigned int segments, double len) {
	Curve c;
	double x = 0.0;
	std::normal_distribution length(len, len / 2);
	double y = length(g);
	while (y <= 0)
		y = length(g);
	assert(y > 0.0);
	double unit = M_PI / 2 / segments;
	double mean = (segments - 0.5) * unit;
	std::normal_distribution angle(mean, mean / 3);
	double r = angle(g); // until r > 0
	while (r <= 0.0 || r >= M_PI / 2)
		r = angle(g);
	assert(r > 0.0);
	assert(r < M_PI / 2.0);
	c += Point{x, y, tan(r)};

	while (--segments > 0) {
		double l = length(g);
		while (l <= 0.0)
			l = length(g);
		assert(l > 0.0);
		x += l*cos(r);
		y += l*sin(r);

		double mean = (segments - 0.5) * unit;
		std::normal_distribution angle(mean, mean / 3);
		double nr = angle(g); // until nr > 0 && nr < r
		while (nr <= 0.0 || nr >= r)
			nr = angle(g);
		assert(nr > 0.0);
		assert(nr < r);
		c.emplace_back(x, y, tan(nr));
		r = nr;
	}
	return c;
}

template<typename... C>
static Point
recalculate(C const&... curves) {
	Point p(0.0, 0.0, 0.0);
	// x    = max x_i
	// y    = sum y_i + r_i * (max x_j - x_i)
	// rate = sum rate_i
	using expander = int[];
	expander{
		(void(
			[&p](const C& curve) {
				const Point& c = curve.back();
				if (c.x > p.x)
					p.x = c.x;
				p.y += c.y - c.rate * c.x;
				p.rate += c.rate;
			}(curves)), 0)...
	};
	p.y += p.rate * p.x;
	return p;
}
