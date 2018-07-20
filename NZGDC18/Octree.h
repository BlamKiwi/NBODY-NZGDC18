#pragma once

#include "stdafx.h"

#include "Vec4.h"
#include <array>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>

namespace brandonpelfrey {

	/**!
	 *
	 */
	class Octree {
		// Physical position/mass.
		Vec4 origin;         //! The physical center of this node
	

		// The tree has up to eight children and can additionally store
		// a point, though in many applications only, the leaves will store data.
		std::array<Octree*, 8> children; //! Pointers to child octants
		bool is_clean;

		/*
				Children follow a predictable pattern to make accesses simple.
				Here, - means less than 'origin' in that dimension, + means greater than.
				child:	0 1 2 3 4 5 6 7
				x:      - - - - + + + +
				y:      - - + + - - + +
				z:      - + - + - + - +
		 */

		public:
		Octree() 
			: origin(Vec4(0.0, 0.0, 0.0, 0.0))
		, is_clean(true){
				// Initially, there are no children
				for(int i=0; i<8; ++i) 
					children[i] = nullptr;
			}

		Octree(Octree&& other)
			: origin(other.origin), children(other.children) {
			for (int i = 0; i < 8; ++i)
				children[i] = nullptr;
			}

		~Octree() {
			// Recursively destroy octants
			for(int i=0; i<8; ++i) 
				delete children[i];
		}

		// Determine which octant of the tree would contain 'point'
		int getOctantContainingPoint(const Vec4& point) const {
			int oct = 0;
			if(point.x >= origin.x) oct |= 4;
			if(point.y >= origin.y) oct |= 2;
			if(point.z >= origin.z) oct |= 1;
			return oct;
		}

		bool isLeafNode() const {

			// We are a leaf if we have no children. Since we either have none, or 
			// all eight, it is sufficient to just check the first.
			return children[0] == nullptr;
		}

		void insert(Vec4 point) {
			// If this node doesn't have a data point yet assigned 
			// and it is a leaf, then we're done!
			if(isLeafNode()) {
				
				// Are we the same point in space?
				if (is_clean) {
					origin = point;
					is_clean = false;
				} else if (point.x == origin.x && point.y == origin.y && point.z == origin.z)
				{
					// Accumulate the masses
					origin.w += point.w;
				}
				
				else {
					// We're at a leaf, but there's already something here
					// We will split this node so that it has 8 child octants
					// and then insert the old data that was here, along with 
					// this new data point
					
					// Split the current node and create new empty trees for each
					// child octant.
					for(int i=0; i<8; ++i) {
						children[i] = new Octree();
					}

					// Calculate new centre of mass
					const Vec4 old = origin;
					origin = CentreofMass(old, point);


					// Re-insert the old point, and insert this new point
					// (We wouldn't need to insert from the root, because we already
					// know it's guaranteed to be in this section of the tree)
					auto oct_origin = getOctantContainingPoint(old);
					auto oct_point = getOctantContainingPoint(point);
					assert(oct_point != oct_origin);
					children[oct_origin]->insert(old);
					children[oct_point]->insert(point);
					UpdateCentreOfMass();
				}
			} else {
				// We are at an interior node. Insert recursively into the 
				// appropriate child octant
				int octant = getOctantContainingPoint(point);
				children[octant]->insert(point);
				UpdateCentreOfMass();
			}
		}

		void getPointsInsideRadiusLeafImpl(const Vec4 & source, double radius_sqr, std::vector<Vec4> & results) const
		{
			const Vec4 diff = source - origin;
			const double dist = diff.normSquared();
			if (dist <= radius_sqr)
			{
				results.push_back(origin);
			}
		}

		void getPointsInsideRaduisInnerImpl(const Vec4 & source, double radius_sqr, std::vector<Vec4> & results) const
		{
			// We're at an interior node of the tree. We will check to see if
			// the query radius lies outside the octants of this node.
			for (int i = 0; i<8; ++i) {
				// Is the centre of mass within the radius of influence
				const Vec4 diff = source - origin;
				const double dist = diff.normSquared();
				if (dist > radius_sqr)
				{
					// Centre of mass is outside influence. Use approximation for cluster. 
					results.push_back(origin);
				}
				else
				{
					children[i]->getPointsInsideRadiusSqr(source, radius_sqr, results);
				}

			}
		}

		void getPointsInsideRadiusSqr(const Vec4& source, double radius_sqr, std::vector<Vec4>& results) const
		{
			// If we're at a leaf node, just see if the current data point is inside
			// the query bounding box
			if (isLeafNode() && !is_clean) {
			
				getPointsInsideRadiusLeafImpl(source, radius_sqr, results);
				
			}
			else {
				getPointsInsideRaduisInnerImpl(source, radius_sqr, results);
			}
		}

		protected:
			static Vec4 CentreofMass(Vec4 a, Vec4 b)
			{
				double x_acc = 0.0;
				double y_acc = 0.0;
				double z_acc = 0.0;
				double w_acc = 0.0;

				x_acc += a.x * a.w;
				y_acc += a.y * a.w;
				z_acc += a.z * a.w;
				w_acc += a.w;

				x_acc += b.x * b.w;
				y_acc += b.y * b.w;
				z_acc += b.z * b.w;
				w_acc += b.w;

				return Vec4(x_acc / w_acc,
				 y_acc / w_acc,
				 z_acc / w_acc,
				 w_acc);
			}

			void UpdateCentreOfMass()
			{
				// Centre of mass can be calculated by the 
				// sum of mass-position products over the total mass of the system
				assert(!isLeafNode());
				double x_acc = 0.0;
				double y_acc = 0.0;
				double z_acc = 0.0;
				double w_acc = 0.0;

				for (auto &c : children)
				{
					const Vec4 p = c->origin;
					x_acc += p.x * p.w;
					y_acc += p.y * p.w;
					z_acc += p.z * p.w;
					w_acc += p.w;
				}

				origin.x = x_acc / w_acc;
				origin.y = y_acc / w_acc;
				origin.z = z_acc / w_acc;
				origin.w = w_acc;
			}
	};

}
