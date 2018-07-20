// NZGDC18.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <chrono>
#include <random>
#include <vector>

#include "Octree.h"
#include "Vec4.h"
#include <chrono>
#include <iostream>

const constexpr size_t POINTS = 1000000;
const constexpr size_t ITERATIONS = 7;
const constexpr double DT = 1.0 / 60.0;
const constexpr double G = 6.67408e-11;
const constexpr double TAU = 0.25;

std::vector<typename Vec4> GeneratePoints()
{
	std::mt19937_64 rand;
	std::uniform_real_distribution<typename Vec4::NumericalT> dist(0.0, 1.0);
	std::vector<Vec4> res;
	for (size_t i = 0; i < POINTS; i++)
	{
		res.push_back(Vec4(dist(rand), dist(rand), dist(rand), dist(rand)));
	}
	return res;
}

brandonpelfrey::Octree ConstructOctTree(const std::vector<typename Vec4>& points)
{
	brandonpelfrey::Octree res;
	for (auto& p : points)
	{
		res.insert(p);
	}
	return res;
}



Vec4 Force(const Vec4& a, const Vec4& b, const double G)
{
	Vec4 offs = b - a;
	double r2 = offs.normSquared();
	Vec4 n = offs.normalized();
	double magnitude = G * ((a.w * b.w) / r2);
	Vec4 force = magnitude * n;
	return force;
}

void Integrate(std::vector<Vec4> frame, const brandonpelfrey::Octree& tree, const double dt, const double G)
{
	std::vector<Vec4> scratch;
	scratch.reserve(POINTS);
	size_t i = 0; 
	for (auto &p : frame)
	{
		Vec4 force(0.0, 0.0, 0.0, 0.0);
		scratch.clear( );
		tree.getPointsInsideRadiusSqr(p, TAU * TAU, scratch);
		for (auto& q : scratch)
		{
			force += Force(p, q, G);
		}
		p += dt * force;
	}
}

int main()
{
	auto frame = GeneratePoints();

	auto p1 = std::chrono::steady_clock::now();

	for (size_t i = 0; i < ITERATIONS; i++)
	{
		auto tree = ConstructOctTree(frame);
		Integrate(frame, tree, DT, G);
	}

	auto p2 = std::chrono::steady_clock::now();

	auto diff = p2 - p1;
	double total = std::chrono::duration_cast<std::chrono::duration<double>>(diff).count();
	double frame_time = total / double(ITERATIONS);
	double fps = 1.0 / frame_time;

	std::cerr << "Average rate for " << POINTS << " points is " << fps << " fps." << std::endl;

    return 0;
}

