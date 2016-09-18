#include "stdafx.h"
#include <algorithm>

inline float Vector2f::length() const
{
	return sqrt(x * x + y * y);
}

inline float Vector2f::operator*(Vector2f v) const
{
	return x * v.x + y * v.y;
}

inline Vector2f Vector2f::operator+(Vector2f v) const
{
	return Vector2f(x + v.x, y + v.y);
}

Vector2f Vector2f::operator-(Vector2f v) const
{
	return Vector2f(x - v.x, y - v.y);
}

inline Vector2f Vector2f::operator*(float f) const
{
	return Vector2f(x * f, y * f);
}

Vector2f Vector2f::normalize() const
{
	Vector2f v = *(this);
	
	float l = v.length();

	v.x /= l;
	v.y /= l;

	return v;
}

void FlowSolver::splitCurve()
{
	if (M < 4)
		throw "Not enough singularities";

	singularities.clear();
	colocations.clear();
	normals.clear();

	/*
		(0.25, 0.5)
		(0.0, 0.5)
		(0.0, 0.0)
		(0.25, -0.25)
		(0.0, -0.5)

		length = 0.25 + 0.5 + Pi * 0.25
	*/

	float l1 = 4.0f/9.0f, l2 = 2.0f/9.0f, l3 = 2.0f/9.0f, l4 = 1.0f/9.0f;
	int m1, m2, m3, m4;
	float L = l1 + l2 + l3 + l4;
	m1 = std::max(1, int(l1 * (M - 1) / L));
	m2 = std::max(1, int(l2 * (M - 1) / L));
	m3 = std::max(1, int(l3 * (M - 1) / L));
	m4 = (M - 1) - m1 - m2 - m3;

	big_delta = std::min(std::min(std::min(l1 / m1, l2 / m2), l3 / m3), l4 / m4);
	delta = big_delta * 0.45f;

	singularities.push_back(Vector2f(0.0f, -2.0f/9.0f));

	for (int i = 0; i < m1; i++)
	{
		float x, y;

		x = 0.0f;
		y = -2.0f/9.0f + (i + 1) * (4.0/9.0f) / m1;
		singularities.push_back(Vector2f(x, y));

		x = 0.0f;
		y = -2.0f/9.0f + (i + 0.5f) * (4.0/9.0f) / m1;
		colocations.push_back(Vector2f(x, y));
		normals.push_back(Vector2f(-1.0f, 0.0f).normalize());
	}

	for (int i = 0; i < m2; i++)
	{
		float x, y;

		x = 0.0f + (i + 1) * (2.0f/9.0f) / m2;
		y = 2.0f/9.0f;
		singularities.push_back(Vector2f(x, y));

		x = 0.0f + (i + 0.5f) * (2.0f/9.0f) / m2;
		y = 2.0f/9.0f;
		colocations.push_back(Vector2f(x, y));
		normals.push_back(Vector2f(0.0f, 1.0f).normalize());
	}

	for (int i = 0; i < m3; i++)
	{
		float x, y;

		x = l2;
		y = l2 - (i + 1) * l2 / m3;
		singularities.push_back(Vector2f(x, y));

		x = l2;
		y = l2 - (i + 0.5f) * l2 / m3;
		colocations.push_back(Vector2f(x, y));
		normals.push_back(Vector2f(1.0f, 0.0f).normalize());
		//normals.push_back(Vector2f(-1.0f, 0.0f).normalize());
	}

	for (int i = 0; i < m4; i++)
	{
		float x, y;

		x = l2 - (i + 1) * l4 / m4;
		y = 0.0f;
		singularities.push_back(Vector2f(x, y));

		x = l2 - (i + 0.5f) * l4 / m4;
		y = 0.0f;
		colocations.push_back(Vector2f(x, y));
		normals.push_back(Vector2f(0.0f, 1.0f).normalize());
		//normals.push_back(Vector2f(1.0f, 0.0f).normalize());
	}

	stage = 1;
	return;
}

float FlowSolver::u_j(int j, float x, float y) const
{
	if (stage < 1)
		throw "Curve is not splitted";
	if (j >= M)
		throw "Point index (j) out of range";

	float x_j = singularities[j].x;
	float y_j = singularities[j].y;
	
	float dist_sq = (x - x_j) * (x - x_j) + (y - y_j) * (y - y_j);
	dist_sq = std::max(dist_sq, delta * delta);

	float res = 1.0f / (2.0f * Pi) * (y_j - y) / dist_sq;

	return res;
}

float FlowSolver::v_j(int j, float x, float y) const
{
	if (stage < 1)
		throw "Curve is not splitted";
	if (j >= M)
		throw "Point index (j) out of range";

	float x_j = singularities[j].x;
	float y_j = singularities[j].y;
	
	float dist_sq = (x - x_j) * (x - x_j) + (y - y_j) * (y - y_j);
	dist_sq = std::max(dist_sq, delta * delta);

	float res = 1.0f / (2.0f * Pi) * (x - x_j) / dist_sq;

	return res;
}

Vector2f FlowSolver::V_j(int j, float x, float y) const
{
	if (stage < 1)
		throw "Curve is not splitted";
	if (j >= M)
		throw "Point index (j) out of range";

	return Vector2f(u_j(j, x, y), v_j(j, x, y));
}

float FlowSolver::phi_j(int j, float x, float y) const
{
	if (stage < 1)
		throw "Curve is not splitted";
	if (j >= M)
		throw "Point index (j) out of range";

	float x_j = singularities[j].x;
	float y_j = singularities[j].y;
	float dist = sqrt((x - x_j) * (x - x_j) + (y - y_j) * (y - y_j));

	if (dist < delta)
		return 0.0f;

	return std::atan2(y - y_j, x - x_j);
}

inline float FlowSolver::u_Inf() const
{
	return cos(alpha);
}

inline float FlowSolver::v_Inf() const
{
	return sin(alpha);
}

inline Vector2f FlowSolver::V_Inf() const
{
	return Vector2f(cos(alpha), sin(alpha));
}

Vector2f FlowSolver::V(float x, float y) const
{
	if (stage < 3)
		throw "System of equations is not solved";
	
	Vector2f res = V_Inf();
	for (int j = 0; j < M; j++)
		res = res + V_j(j, x, y) * gamma[j];

	return res;
}

float FlowSolver::phi(float x, float y) const
{
	/* TODO: modify this function to avoid stripped output*/
	if (stage < 3)
		throw "System of equations is not solved";

	// // VORTEX METHODE 
	//float res = x * u_Inf() + y * v_Inf();
	//for (int j = 0; j < M; j++)
	//	res += gamma[j] * phi_j(j, x, y)/* / 2.0f / Pi*/;

	float res = 0.0f;
	
	float gamma_sum = 0.0f;
	compf tmp_res(0.0f, 0.0f);
	compf z(x, y);
	
	for (int j = 0; j < M - 1; j++)
	{
		gamma_sum += gamma[j];
		compf z01(singularities[j].x, singularities[j].y);
		compf z02(singularities[j + 1].x, singularities[j + 1].y);
		compf z_star = (z01 + z02) / 2.0f;
		compf A = gamma_sum * (z02 - z01);

		if (abs(z - z_star) > delta)
			tmp_res += A / (2.0f * Pi * compf(0.0, 1.0) * (z - z_star));
	}
	
	res = tmp_res.real();
	res += gamma0 * phi_j(M - 1, x, y) / 2.0f / Pi;
	res += x * u_Inf() + y * v_Inf();

	return res;
}

float FlowSolver::psi(float x, float y) const
{
	if (stage < 3)
		throw "System of equations is not solved";
	
	float res = exp(y * u_Inf() - x * v_Inf());;
	for (int j = 0; j < M; j++)
	{
		float x_j = singularities[j].x;
		float y_j = singularities[j].y;
	
		float dist = sqrt((x - x_j) * (x - x_j) + (y - y_j) * (y - y_j));
		dist = std::max(dist, delta);

		res *= pow(dist, - gamma[j] / (2.0f * Pi));
	}

	res = log(res);

	return res;
}

float FlowSolver::Cp(float x, float y) const
{
	if (stage < 3)
		throw "System of equations is not solved";

	Vector2f v1 = V(x, y);
	Vector2f v2 = V_Inf();

	float v1_sq = v1.x * v1.x + v1.y * v1.y;
	float v2_sq = v2.x * v2.x + v2.y * v2.y;
	float res = 1.0f - v1_sq / v2_sq;

	return res;
}

void FlowSolver::setGamma0(float value)
{
	gamma0 = value;
	stage = std::min(stage, 1);
}

void FlowSolver::setAlpha(float value)
{
	if (!(- Pi / 2.0f < alpha && alpha < Pi / 2.0f))
		throw "New value of alpha must be from (-Pi/2; Pi/2)";

	alpha = value;
	stage = std::min(stage, 1);
}

void FlowSolver::setM(int value)
{
	if (!(4 <= value && value <= 1000))
		throw "New value of M must be from [4; 100]";

	M = value;
	stage = std::min(stage, 0);
}

void FlowSolver::makeGammaSystem()
{
	if (stage < 1)
		throw "Curve is not splitted";

	A.clear();
	A.resize(M, std::vector<float>(M + 1));
	
	for (int k = 0; k < M - 1; k++)
	{
		for (int j = 0; j < M; j++)
		{
			A[k][j] = V_j(j, colocations[k].x, colocations[k].y) * normals[k];
		}
		A[k][M] = - (V_Inf() * normals[k]);
	}

	for (int j = 0; j < M; j++)
	{
		A[M - 1][j] = 1.0f;
	}
	A[M - 1][M] = gamma0;

	stage = 2;
	return;
}

void FlowSolver::solveGamma()
{
	if (stage < 2)
		throw "System of equations is not created";

	gamma.clear();

	for (int k = 0; k < M; k++)
	{
		double mx = -1;
		int v = -1;
		for (int l = k; l < M; l++)
		{
			if (abs(A[l][k]) > mx) {
				mx = abs(A[l][k]);
				v = l;
			}
		}
		if (v != k) swap(A[v], A[k]);

		if (abs(A[k][k]) <= Eps)
			throw "System of equations is singular";

		for (int j = k + 1; j <= M; j++)
			A[k][j] /= A[k][k];
		A[k][k] = 1.0f;

		for (int l = 0; l < M; l++)
			if (l != k)
			{
				for (int j = k + 1; j <= M; j++)
					A[l][j] -= A[k][j] * A[l][k];
				A[l][k] = 0.0f;
			}
	}

	for (int k = 0; k < M; k++)
		gamma.push_back(A[k][M]);

	stage = 3;
	return;
}

void FlowSolver::solve()
{
	splitCurve();
	makeGammaSystem();
	solveGamma();
}

std::vector<Vector2f> FlowSolver::getColocations() const
{
	return colocations;
}

std::vector<Vector2f> FlowSolver::getSingularities() const
{
	return singularities;
}

float FlowSolver::getDelta() const
{
	return delta;
}

float FlowSolver::V_scal(float x, float y) const
{
	Vector2f v = V(x, y);
	return sqrt(v.x * v.x + v.y * v.y);
}
