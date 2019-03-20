#include "pch.h"
#include <iostream>
#include <fstream>
#include <cstdio>

using namespace std;

void height_analyz(int kratnost, double MASS_Fa_nor[28][2], double m, double g, long Ho, int f, double *TIMES) {

	double Vmax = 0, Fa_max, Hn, H_buf = 0, TIME = 0, time_buf, V0 = 0;
	int V = 10, n, k, s = 0, p = 0;
	long Hd = 0;
	double *MASS_Fa = new double[f];

	for (int i = 0; i < f; i++) {
		MASS_Fa[i] = MASS_Fa_nor[i][1];
	}

	Fa_max = g * m;
	//Рассчет максимальной скорости
	for (int i = 0; i < 28; i++) {
		if (MASS_Fa[i] <= Fa_max && Fa_max <= MASS_Fa[i + 1]) {
			if ((abs(MASS_Fa[i] - Fa_max)) >= (abs(MASS_Fa[i + 1] - Fa_max))) {
				Vmax = V * (i + 1 - abs((Fa_max - MASS_Fa[i + 1]) / (MASS_Fa[i + 1] - MASS_Fa[i])));
			}
			else {
				Vmax = V * (i + 1 - abs((Fa_max - MASS_Fa[i]) / (MASS_Fa[i + 1] - MASS_Fa[i])));
			}
			n = i + 1;
			break;
		}
	}
	for (int i = 0; i < f - 1; i++) {
		MASS_Fa[i] = (MASS_Fa[i] + MASS_Fa[i + 1]) / 2;
	}
	for (int i = 0; i < f; i++) {
		MASS_Fa[i] /= m;
	}
	//Рассчет высоты для набора максимальной скорости
	double *height = new double[n];
	for (int i = 1; i < n; i++) {
		height[i - 1] = V * V*(2 * i - 1) / 2 / (g - MASS_Fa[i - 1]);
		Hd += height[i - 1];
	}
	height[n - 1] = (Vmax*Vmax - V * V*(n - 1)*(n - 1)) / 2 / (g - MASS_Fa[n - 1]);
	Hd += height[n - 1];
	//Анализ высоты и построение таблицы времени
	k = Hd / kratnost;
	double *times = new double[n]; //время по перемене сопротивлений
	for (int i = 0; i < k; i++) {
		if (i == 0) {
			while (H_buf < kratnost) {
				s++;
				H_buf += height[p];
				p++;
			}
			s--;
			p--;
			H_buf -= height[p];
			times[0] = V / (g - MASS_Fa[0]);
			for (int i = 1; i < s + 1; i++) {
				times[i] = V / (g - MASS_Fa[i]);
			}
			for (int i = 0; i < s; i++) {
				TIME += times[i];
			}
			time_buf = (sqrt(2 * (kratnost*(i + 1) - H_buf)*(g - MASS_Fa[s]) + V * V*s*s) - s * V) / (g - MASS_Fa[s]);
			TIME += time_buf;
			TIMES[i] = TIME;
			time_buf = times[s] - time_buf;
			s = 0;
			if ((i + 1)*kratnost == Ho) { break; }
		}
		else {
			if (kratnost < (V*(i + 1) + (g - MASS_Fa[p])*time_buf / 2)*time_buf) {
				V0 = Vmax - time_buf * (g - MASS_Fa[p]);
				Hn = sqrt(2 * kratnost * (g - MASS_Fa[p]) + V0 * V0);
				TIMES[i] = (Hn - V0) / (g - MASS_Fa[p]);
				if ((i + 1)*kratnost == Ho) { break; }
				i++;
				V0 = Hn;
				Hn = sqrt(2 * kratnost * (g - MASS_Fa[p]) + V0 * V0);
				TIMES[i] = (Vmax - Hn) / (g - MASS_Fa[p]);
				if ((i + 1)*kratnost == Ho) { break; }
				H_buf = (Vmax*Vmax - Hn * Hn) / 2 / (g - MASS_Fa[p]);
				TIMES[i] += (kratnost - H_buf) / Vmax;
				if ((i + 1)*kratnost == Ho) { break; }
				time_buf = kratnost / Vmax;
				while (time_buf != -1) {
					i++;
					TIMES[i] = time_buf;
					if ((i + 1)*kratnost == Ho) { break; }
				}
			}
			else {
				TIMES[i] = time_buf;
				while (H_buf < (i + 1)*kratnost) {
					s++;
					H_buf += height[p];
					p++;
				}
				s--;
				p--;
				H_buf -= height[p];
				times[p] = V / (g - MASS_Fa[p]);
				for (int i = p; i < p + s + 1; i++) {
					times[i] = V / (g - MASS_Fa[i]);
				}
				time_buf = (sqrt(2 * (kratnost*(i + 1) - H_buf)*(g - MASS_Fa[p]) + V * V*(p)*(p)) - p * V) / (g - MASS_Fa[p]);
				TIMES[i] += time_buf;
				time_buf = times[p] - time_buf; s = 0;
				if ((i + 1)*kratnost == Ho) { break; }
			}
		}
	}
}

void func(const double visota, double veter[15][3], int kvet, double *time, const double Vo, double Fa[28][2], int kFa, double m, double th, double &a, double &vs, double &z) {
	int i = 0, j = 0, kratnost = Fa[1][0] - Fa[0][0];
	double Wx = 0, x = 0, Wz = 0, to = 0, sum1 = 0, sum2 = 0, temp = 0, Fsr = 0, ayck, x_f=0;
	bool tm = 0;
	int d = 0;
	int act = 0;
	double copytemp = 0;
	a = 0;
	vs = 0;
	z = 0;
	for (i = kvet - 1; visota < veter[i][0]; i--) {

	};
	double **tabl = new double*[25];
	for (j = 0; j < 25; j++) {
		tabl[j] = new double[4];
	}
	int tempi = i;
	for (i; i > 0; i--) {
		Wx = (veter[i][1] + veter[i - 1][1]) / 2;
		x = x + Wx * time[i - 1];
		Wz = (veter[i][2] + veter[i - 1][2]) / 2;
		z = z + Wz * time[i - 1];
	}
	for (j = kFa - 1; Vo < Fa[j][0]; j--) {

	};
	int tempj = j;
	int sosi2 = 0;
	for (j; j > 0; j--) {
		Fsr = (Fa[j][1] + Fa[j - 1][1]) / 2;
		temp = kratnost * m / Fsr;
		to = to + temp;
		if (to > th) {
			tm = 1;
			to = temp - to + th;
			j++;
			break;
		}
		sum1 = sum1 + 1 / Fsr;
		sum2 = sum2 + sosi2 / Fsr;
		sosi2++;
	}
	if (tm == 1) {
		vs = Fa[j][0] * to - ((Fa[j - 1][1] - ((to*(Fa[j - 1][1] - Fa[j][1])) / (2 * temp))) * to * to) / (2 * m);
	}
	vs = vs + kratnost * m*((Vo - kratnost / 2)*sum1 - kratnost * sum2);
	vs = vs + x;
	a = x / (sqrt(x*x + z * z));

	a = acos(a);
	cout << endl << "Угол: " << a * 57.2958 << endl;
	vs = z * sin(a) + vs * cos(a);
	z = z * cos(a) - vs * sin(a);
	x_f = vs;
	cout << "Координаты (x, z): " << vs << ", " << z << endl;
	to = 0;
	////////////////////////////////////////////////////////////////////////////////
	d = 0;
	j = tempj;
	int p = j;
	int sosi = 0;
		int speed = Vo;
	ofstream coord("Coords.txt");
	coord << "z;x;V;t" << endl;
	coord << z <<";"<< x_f <<";"<< Vo<<";" << "0";
	for (i = j; i > 0; i--) {
		Fsr = (Fa[i][1] + Fa[i - 1][1]) / 2;
		temp = kratnost * m / Fsr;
		ayck = Fsr / m;
		to += temp;
		if (to > th) {
			to -= temp;
			to = th - to;
			ayck = (Fa[p - 1][1] - to * (Fa[p - 1][1] - Fa[p][1]) / 2 / temp) / m;
			vs = Fa[p][0] * to - ayck*to*to / 2;
			tabl[d][0] = 0;
			tabl[d][1] = 0;
			tabl[d][2] = speed - ayck*to;
			coord << endl << tabl[d][0] << ";" << tabl[d][1] << ";" << tabl[d][2] << ";" << th;
			break;
		}
		vs = 0;
		vs = (Fa[p][0] * Fa[p][0] - Fa[p - 1][0] * Fa[p - 1][0]) / ayck / 2;



		if (d > 0) {
			z = z + sin(a)*vs;
			x_f = x_f  - cos(a)*vs;
			tabl[d][0] = z;
			tabl[d][1] = x_f;
		}
		else {
			z = z + sin(a)*vs;
			x_f = x_f - cos(a)*vs;
			tabl[d][0] = z;
			tabl[d][1] = x_f;
		}

		act = 0;
		copytemp = temp;
		do {
			if (copytemp < time[tempi - 1]) {
				Wx = (veter[tempi][1] + veter[tempi - 1][1]) / 2;
				tabl[d][1] = tabl[d][1] + Wx * copytemp;
				Wz = (veter[tempi][2] + veter[tempi - 1][2]) / 2;
				tabl[d][0] = tabl[d][0] + Wz * copytemp;
				time[tempi - 1] = time[tempi - 1] - copytemp;
				act = 1;
			}
			else {
				if (copytemp >= time[tempi - 1]) {
					Wx = (veter[tempi][1] + veter[tempi - 1][1]) / 2;
					tabl[d][1] = tabl[d][1] + Wx * time[tempi - 1];
					Wz = (veter[tempi][2] + veter[tempi - 1][2]) / 2;
					tabl[d][0] = tabl[d][0] + Wz * time[tempi - 1];
					tempi--;
					copytemp = copytemp - time[tempi];
				}
			}
		} while (act == 0);

		tabl[d][2] = speed - 10;
		tempj = j;
		coord << endl << tabl[d][0] << ";" << tabl[d][1] << ";" << tabl[d][2] << ";" << to;
		d++;
		p--;
		speed -= 10;
	}
	coord.close();
}

void main() {

	double m, g = 9.81, MASS_Fa[28][2], veter[15][3], V0 = 0, a = 0, x = 0, z = 0;
	int kratnost = 100, f = 0;
	long Ho;
	char a1;
	setlocale(LC_ALL, "Russian");

	cout << "Введите положительную высоту Ho: ";
	cin >> Ho;
	while (Ho < 100 && Ho>1400) {
		cout << "Значение высоты некорректно, повторите ввод: ";
		cin >> Ho;
	}
	cout << "Введите массу груза: ";
	cin >> m;
	while (m <= 0) {
		cout << "Значение массы некорректно, повторите ввод: ";
		cin >> m;
	}
	cout << "Введите начальную скорость: ";
	cin >> V0;
	while (V0 < 10 && V0>250) {
		cout << "Значение скорости некорректно, повторите ввод: ";
		cin >> V0;
	}

	ifstream Fa("F.csv");
	Fa.imbue(locale("rus_rus.1251"));
	for (int i = 0, j = 0; i < 28; i++) {
		Fa >> MASS_Fa[i][j];
		Fa >> a1;
		j++;
		Fa >> MASS_Fa[i][j];
		j--;
	}
	Fa.close();

	ifstream vet("Wind.csv");
	vet.imbue(locale("rus_rus.1251"));
	for (int i = 0, j = 0; i < 15; i++) {
		vet >> veter[i][j];
		vet >> a1;
		j++;
		vet >> veter[i][j];
		vet >> a1;
		j++;
		vet >> veter[i][j];
		j -= 2;
	}
	vet.close();

	f = Ho / kratnost;
	double *TIMES = new double[f];

	cout << "Введите точное g: ";
	cin >> g;
	while (g < 8 && g>11) {
		cout << "Значение g некорректно, повторите ввод: ";
		cin >> g;
	}
	height_analyz(kratnost, MASS_Fa, m, g, Ho, f, TIMES);

	for (int i = 0; i < f; i++) {
		g += TIMES[i];
	}

	double *TIMES_INVERT = new double[f];

	for (int i = 0; i < f; i++) {
		TIMES_INVERT[i] = TIMES[f - i - 1];
	}

	func(Ho, veter, 15, TIMES_INVERT, V0, MASS_Fa, 28, m, g, a, x, z);

	cout << "Значения траектории выгружены в файл Coords.txt. Траектория является примерной из-за чего может наблюдаться несоответствие координатам выше." << endl;

	system("pause");

	delete[] TIMES;
	delete[] TIMES_INVERT;
}
