#pragma once

#include <vector>
#include <cmath>
#include <complex>

typedef std::complex<float> compf;

const float Pi = 3.1415926535f;
const float Eps = 1e-7f;
const float Inf = 1e7f;

struct Vector2f
{
	float x, y;

	Vector2f(float x = 0.0f, float y = 0.0f): x(x), y(y) {}

	inline float length() const;
	inline float operator*(Vector2f v) const;
	inline Vector2f operator+(Vector2f v) const;
	Vector2f operator-(Vector2f v) const; // TODO: why cannot be implemented as inline?
	inline Vector2f operator*(float f) const;
	Vector2f normalize() const;
};

struct FlowSolver
{
private:
	float gamma0;
	float alpha;
	int M;
	float delta, big_delta;

	int stage;
	std::vector<Vector2f> colocations;
	std::vector<Vector2f> singularities;
	std::vector<Vector2f> normals;
	std::vector< std::vector<float> > A;
	std::vector< float > gamma;

	float u_j(int j, float x, float y) const;
	float v_j(int j, float x, float y) const;
	Vector2f V_j(int j, float x, float y) const;
	inline float u_Inf() const;
	inline float v_Inf() const;
	inline Vector2f V_Inf() const;
	float phi_j(int j, float x, float y) const;

public:
	FlowSolver(): stage(0), gamma0(1.0f), alpha(0.0f), M(10) {};
	void setGamma0(float gamma0);
	void setAlpha(float alpha);
	void setM(int M);
	void splitCurve();
	void makeGammaSystem();
	void solveGamma();
	void solve();

	std::vector<Vector2f> getColocations() const;
	std::vector<Vector2f> getSingularities() const;
	float getDelta() const;

	Vector2f V(float x, float y) const;
	float Cp(float x, float y) const;
	float phi(float x, float y) const;
	float psi(float x, float y) const;
	float V_scal(float x, float y) const;
};
