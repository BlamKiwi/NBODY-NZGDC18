#pragma once

#include <cmath>

template <typename F>
struct Vector4;

template <typename F>
Vector4<F> operator*(F r, const Vector4<F>& v);

using Vec4 = Vector4<float>;

/*
	3D Position + Mass
*/
template <typename F>
struct alignas(sizeof(F) * 4) Vector4 {
	using NumericalT = F;

	union {
		struct {
			F x;
			F y;
			F z; 
			F w;
		};
		F D[4];
	};

	Vector4() { }
	Vector4(F _x, F _y, F _z, F _w)
		:x(_x), y(_y), z(_z), w(_w)
	{ }

	F& operator[](unsigned int i) {
		return D[i];
	}

	const F& operator[](unsigned int i) const {
		return D[i];
	}

	F maxComponent() const {
		F r = x;
		if(y>r) r = y;
		if(z>r) r = z;
		return r;
	}

	F minComponent() const {
		F r = x;
		if(y<r) r = y;
		if(z<r) r = z;
		return r;
	}

	Vector4 operator+(const Vector4& r) const {
		return Vector4(x+r.x, y+r.y, z+r.z, F()); 
	}

	Vector4 operator-(const Vector4& r) const {
		return Vector4(x-r.x, y-r.y, z-r.z, F()); 
	}

	Vector4 cmul(const Vector4& r) const {
		return Vector4(x*r.x, y*r.y, z*r.z, F());
	}

	Vector4 cdiv(const Vector4& r) const {
		return Vector4(x/r.x, y/r.y, z/r.z, F());
	}

	Vector4 operator*(F r) const {
		return Vector4(x*r,y*r,z*r, F());
	}


	Vector4 operator/(F r) const {
		return Vector4(x/r, y/r, z/r, F());
	}

	Vector4& operator+=(const Vector4& r) {
		x+=r.x;
		y+=r.y;
		z+=r.z;
		w+=r.w;
		return *this;
	}

	Vector4& operator-=(const Vector4& r) {
		x-=r.x;
		y-=r.y;
		z-=r.z;
		w-=r.w;
		return *this;
	}

	Vector4& operator*=(F r) {
		x *= r; y *= r; z *= r; w *= r;
		return *this;
	}

	// Inner/dot product
	F operator*(const Vector4& r) const {
		return x*r.x + y*r.y + z*r.z;
	}

	F norm() const {
		return sqrtf(x*x+y*y+z*z);
	}

	F normSquared() const {
		return x*x + y*y + z*z;
	}

	// Cross product
	Vector4 operator^(const Vector4& r) const {
		return Vector4(
				y * r.z - z * r.y, 
				z * r.x - x * r.z, 
				x * r.y - y * r.x,
				F()
				);
	}

	Vector4 normalized() const {
		return *this / norm();
	}
};

template <typename F>
inline Vector4<F> operator*(F r, const Vector4<F>& v) {
	return Vector4<F>(v.x*r, v.y*r, v.z*r, v.w*r);
}

