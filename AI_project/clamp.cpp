#include "clamp.h"


void clamp_func::clamp_sy(double** matrix, int col, int row, double min, double max)
{
	int i, j;
	for (i = 0; i < col; i++)
	{
		for (j = 0; j < row; j++)
		{
			if (matrix[i][j] > max)
			{
				matrix[i][j] = max;
			}
			else if (matrix[i][j] < min)
			{
				matrix[i][j] = min;
			}
		}
	}
}

