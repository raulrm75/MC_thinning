#define ROWS %(rows)d
#define COLS %(cols)d

__device__ int dimension(int r, int c)
{
	return (r %% 2) + (c %% 2);
}
__global__ void sumFacets(int *K, int *cF)
{
	int row = blockIdx.x * blockDim.x + threadIdx.x;
	int col = blockIdx.y * blockDim.x + threadIdx.y;
	int idx = row * COLS + col;
	
	int dim, r,c,d,ix;
	
	if((row < ROWS) && (col < COLS))
	{
		dim = dimension(row, col);
		int facets = 0;
		for(int dr = -1; dr <= 1; dr++)
		{
			for(int dc = -1; dc <= 1; dc++)
			{
				if(!((dr != 0) && (dc != 0)))
				{
					r = row + dr;
					c = col + dc;
					ix = r * COLS + c;
					if((0 <= r) && (r < ROWS) && (0 <= c) && (c < COLS))
					{
						d = dimension(r, c);
						if(((d - dim) == 1) && K[ix] > 0)
						{
							facets++;
						}
					}
				}
			}
		}
		cF[idx] = facets;
	}
}

__global__ void markIsolatedCells(int *K, int *cF, int *mI, int *newIsolated)
{
	int row = blockIdx.x * blockDim.x + threadIdx.x;
	int col = blockIdx.y * blockDim.x + threadIdx.y;
	int idx = row * COLS + col;
	
	if((row < ROWS) && (col < COLS))
	{
		if(cF[idx] == 0)
		{
			if(mI[idx] == 0)
			{
				mI[idx] = 1;
				newIsolated[idx] = 1;
			}
		}					
	}
}

__global__ void getSimpleFacets(int *K, int *cF, int *sFacets)
{
	int row = blockIdx.x * blockDim.x + threadIdx.x;
	int col = blockIdx.y * blockDim.x + threadIdx.y;
	int idx = row * COLS + col;
	
	int dim, r,c,d,ix;
	
	if((row < ROWS) && (col < COLS) && (cF[idx] == 1))
	{
		dim = dimension(row, col);

		for(int i = -1; i <= 1; i++)
		{
			for(int j = -1; j <= 1; j++)
			{
				if(!((i != 0) && (j != 0)))
				{
					r = row + i;
					c = col + j;
					ix = r * COLS + c;
					if((0 <= r) && (r < ROWS) && (0 <= c) && (c < COLS))
					{
						d = dimension(r, c);
						if(((d - dim) == 1) && K[ix] > 0)
						{
							sFacets[ix] = 1;
						}
					}
				}
			}
		}
	}
}

__global__ void getSimpleFaces(int *cF, int * sF, int *sFaces)
{
	int row = blockIdx.x * blockDim.x + threadIdx.x;
	int col = blockIdx.y * blockDim.x + threadIdx.y;
	int idx = row * COLS + col;
	
	int dim, r,c,d,ix, last_idx;
	
	if((row < ROWS) && (col < COLS) && (sF[idx] == 1))
	{
		dim = dimension(row, col);
		last_idx = -1;
		
		for(int i = -1; i <= 1; i++)
		{
			for(int j = -1; j <= 1; j++)
			{
				if(!((i != 0) && (j != 0)))
				{
					r = row + i;
					c = col + j;
					ix = r * COLS + c;
					if((0 <= r) && (r < ROWS) && (0 <= c) && (c < COLS))
					{
						d = dimension(r, c);
						if(((dim - d) == 1) && cF[ix] == 1)
						{
							if(ix > last_idx)
							{
								last_idx = ix;
							}
						}
					}
				}
			}
		}
		sFaces[last_idx] = 1;
	}
}

__global__ void updateR(int *K, int *cR)
{
	int row = blockIdx.x * blockDim.x + threadIdx.x;
	int col = blockIdx.y * blockDim.x + threadIdx.y;
	int idx = row * COLS + col;
	
	if((row < ROWS) && (col < COLS) && (K[idx] == 1))
	{
		cR[idx]++;
	}
}

__global__ void updateI(int *newIsolated, int *cR, int *cI)
{
	int row = blockIdx.x * blockDim.x + threadIdx.x;
	int col = blockIdx.y * blockDim.x + threadIdx.y;
	int idx = row * COLS + col;
	
	if((row < ROWS) && (col < COLS) && (newIsolated[idx] == 1))
	{
		cI[idx] = cR[idx];
	}
}

__global__ void calcMedialPersistence(int *cR, int *cI, int *MPabs, float *MPrel)
{
	int row = blockIdx.x * blockDim.x + threadIdx.x;
	int col = blockIdx.y * blockDim.x + threadIdx.y;
	int idx = row * COLS + col;
	
	if((row < ROWS) && (col < COLS) && (cR[idx] > 0))
	{
		MPabs[idx] = cR[idx] - cI[idx];
		MPrel[idx] = 1. - cI[idx] / cR[idx];
	}
}

__global__ void markRemove(int *sFacets, int *sFaces, int *cI,
	int *MPabs, int *MPrel, int absT, float relT, int useMP, 
	int *rFacets, int *rFaces, int *sR)
{
	int row = blockIdx.x * blockDim.x + threadIdx.x;
	int col = blockIdx.y * blockDim.y + threadIdx.y;
	int idx = row * COLS + col;
	
	if(row < ROWS && col < COLS && sFacets[idx] == 1)
	{
		if(useMP)
		{
			int dim = dimension(row, col);
			if(cI[idx] == 0 || MPabs[idx] < absT || MPrel[idx] < relT)
			{
				rFacets[idx] = 1;
				sR[0] = 1;
			}
		}
		else
		{
			rFacets[idx] = 1;
			sR[0] = 1;
		}
		
		for(int dr = -1; dr <= 1; dr++)
		{
			for(int dc = -1; dc <= 1; dc++)
			{
				if(!(dr == 0 && dc == 0))
				{
					int row2 = row + dr;
					int col2 = col + dc;
					int idx2 = row2 * COLS + col2;
					if(0 <= row2 && row2 < ROWS && 0 <= col2 && col2 < COLS && sFaces[idx2] == 1)
					{
						if(useMP)
						{
							int dim = dimension(row, col);
							if(cI[idx] == 0 || MPabs[idx] < absT || MPrel[idx] < relT)
							{
								rFaces[idx2] = 1;
							}
						}
						else
						{
							rFaces[idx2] = 1;
						}
					}
				}
			}
		}
	}
}

__global__ void doRemove(int *rFacets, int *rFaces, int *K)
{
	int row = blockIdx.x * blockDim.x + threadIdx.x;
	int col = blockIdx.y * blockDim.x + threadIdx.y;
	int idx = row * COLS + col;
	
	if((row < ROWS) && (col < COLS) && (rFacets[idx] > 0 || rFaces[idx] > 0))
	{
		K[idx] = 0;
	}
}

/*
__global__ thinning(int *K, *cR, *cI, *mI, *cF, int absT, float relT, int *halt)
{
	int2 cell, fcell, bcell, maxfcell, maxbcell, d;
	int cIdx, fIdx, cIdx, maxfIdx = 0, maxbIdx = 0;
	int dim, fdim, bdim;
	int isSignificative = 0, MPabs;
	float MPrel;
	
	cell.x = blockIdx.x * blockDim.x + threadIdx.x;
	cell.y = blockIdx.x * blockDim.y + threadIdx.y;
	cIdx  = cell.x * COLS + cell.y;
	
	if(cell.x < ROWS && cell.y < COLS && K[cIdx] > 0)
	{
		// Marcar pares simples 
		// maxfcell es la faceta 
		// maxbcell es la cara 
		if(cF[cIdx] == 1)
		{
			for(d.x = -1; d.x <= 1; d.x++)
			{
				for(d.y = -1; d.y <= 1; d.y++)
				{
					fcell.x = cell.x + d.x;
					fcell.y = cell.y + d.y;
					bIdx = fcell.x * COLS + fcell.y;
					fdim = fcell.x %% 2 + fcell.y %% 2;
					
					if(abs(d.x) + abs(d.y) > 0 && fdim == dim + 1 && K[bIdx] > 0)
					{
						if(fIdx >= maxfIdx)
						{
							maxfcell.x = fcell.x;
							maxfcell.y = fcell.y;
							maxfIdx = fIdx;
						}
					}
				}
			}
			
			// Calcular medidas de persistencia
			MPabs = mR[maxfIdx] - I[maxfIdx];
			MPrel = 1 - I[maxfIdx]/R[maxfIdx];
			fdim = maxfcell.x %% 2 + maxfcell.y %% 2;
			
			isSignificative = ((absT != 0) && (MPabs >= absT)) || ((relT < 1) && (MPrel >= relT));
			
			// Borrar pares simples poco siginificativos
			if(! isSignificative)
			{
				K[maxfIdx] = 0;
				K[maxbIdx] = 0;
				halt[maxfIdx] = 1;
			}
		}
	}
}
*/
