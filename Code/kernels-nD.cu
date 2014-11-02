#define N1 {N1}
#define N2 {N2}
// En las imagenes 2D, N3 = 1
#define N3 {N3}

#define N N1*N2*N3

#define D1 N2*N3
#define D2 N3
#define D3 1

__device__ int index(int3 cell)
{{
	return i.x * N2 * N3 + i.y * N3 + i.z;
}}

__device__ int3 cell(int idx)
{{
	int3 i;
	
	i.x = int(floor(idx / (N2 * N3))) % N1;
	i.y = int(floor(idx / N3)) % N2;
	i.z = idx % N3;
	return i;
}}

__device__ int dimension(int idx)
{{
	// l es la tupla "compilada" y d es la dimension del espacio
	// devuelve la dimension de la tupla "desensamblada"
	int3 i = cell(idx);
	return (i.x % 2) + (i.y % 2) + (i.z % 2);
}}

__global__ void sumFacets(int *K, int *cF)
{{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idx2;
	int dim, r,c,d,ix;
	
	if(idx < N)
	{{
		if(K[idx] > 0)
		{{
			dim = dimension(idx);
			int facets = 0;
			
			idx2 = idx + D1;
			if(idx2 < N && dimension(idx2) == dim + 1 && K[idx2] > 0)
			{{
				facets++;
			}}
			idx2 = idx - D1;
			if(idx2 >= 0 && dimension(idx2) == dim + 1 && K[idx2] > 0)
			{{
				facets++;
			}}
			
			idx2 = idx + D2;
			if(idx2 < N && dimension(idx2) == dim + 1 && K[idx2] > 0)
			{{
				facets++;
			}}
			idx2 = idx - D2;
			if(idx2 >= 0 && dimension(idx2) == dim + 1 && K[idx2] > 0)
			{{
				facets++;
			}}
			
			if(N3 > 1)
			{{
				idx2 = idx + D3;
				if(idx2 < N && dimension(idx2) == dim + 1 && K[idx2] > 0)
				{{
					facets++;
				}}
				idx2 = idx - D3;
				if(idx2 >= 0 && dimension(idx2) == dim + 1 && K[idx2] > 0)
				{{
					facets++;
				}}
			}}
			cF[idx] = facets;
		}}
	}}
}}

__global__ void markIsolatedCells(int *K, int *cF, int *mI, int *newIsolated)
{{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(idx < N && K[idx] > 0)
	{{
		if(cF[idx] == 0)
		{{
			if(mI[idx] == 0)
			{{
				mI[idx] = 1;
				newIsolated[idx] = 1;
			}}
		}}					
	}}
}}

__global__ void getSimpleFacets(int *K, int *cF, int *sFacets)
{{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idx2;
	int dim;
	
	if(idx < N && cF[idx] == 1 && K[idx] > 0)
	{{
		dim = dimension(idx);

		idx2 = idx + D1;
		if(idx2 < N && dimension(idx2) == dim + 1 && K[idx2] > 0)
		{{
			sFacets[idx2] = 1;
		}}
		idx2 = idx - D1;
		if(idx2 >= 0 && dimension(idx2) == dim + 1 && K[idx2] > 0)
		{{
			sFacets[idx2] = 1;
		}}
		
		idx2 = idx + D2;
		if(idx2 < N && dimension(idx2) == dim + 1 && K[idx2] > 0)
		{{
			sFacets[idx2] = 1;
		}}
		idx2 = idx - D2;
		if(idx2 >= 0 && dimension(idx2) == dim + 1 && K[idx2] > 0)
		{{
			sFacets[idx2] = 1;
		}}
		if(N3 > 1)
		{{
			idx2 = idx + D3;
			if(idx2 < N && dimension(idx2) == dim + 1 && K[idx2] > 0)
			{{
				sFacets[idx2] = 1;
			}}
			idx2 = idx - D3;
			if(idx2 >= 0 && dimension(idx2) == dim + 1 && K[idx2] > 0)
			{{
				sFacets[idx2] = 1;
			}}
		}}
	}}
}}

__global__ void getSimpleFaces(int *cF, int * sF, int *sFaces)
{{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	int dim, idx2, dim2, last_idx;
	
	if(idx < N && (sF[idx] == 1))
	{{
		dim = dimension(idx);
		last_idx = -1;
		
		idx2 = idx + D1;
		if(idx2 < N && dimension(idx2) == dim - 1 && cF[idx2] == 1)
		{{
			if(idx2 > last_idx)
			{{
				last_idx = idx2;
			}}
		}}
		idx2 = idx - D1;
		if(idx2 >= 0 && dimension(idx2) == dim - 1 && cF[idx2] == 1)
		{{
			if(idx2 > last_idx)
			{{
				last_idx = idx2;
			}}
		}}
		
		idx2 = idx + D2;
		if(idx2 < N && dimension(idx2) == dim - 1 && cF[idx2] == 1)
		{{
			if(idx2 > last_idx)
			{{
				last_idx = idx2;
			}}
		}}
		idx2 = idx - D2;
		if(idx2 >= 0 && dimension(idx2) == dim - 1 && cF[idx2] == 1)
		{{
			if(idx2 > last_idx)
			{{
				last_idx = idx2;
			}}
		}}
		
		if(N3 > 1)
		{{
			idx2 = idx + D3;
			if(idx2 < N && dimension(idx2) == dim - 1 && cF[idx2] == 1)
			{{
				if(idx2 > last_idx)
				{{
					last_idx = idx2;
				}}
			}}
			idx2 = idx - D3;
			if(idx2 >= 0 && dimension(idx2) == dim - 1 && cF[idx2] == 1)
			{{
				if(idx2 > last_idx)
				{{
					last_idx = idx2;
				}}
			}}
		}}
		sFaces[last_idx] = 1;
	}}
}}

__global__ void updateR(int *K, int *cR)
{{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(idx < N && K[idx] == 1)
	{{
		cR[idx]++;
	}}
}}

__global__ void updateI(int *newIsolated, int *cR, int *cI)
{{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(idx < N && newIsolated[idx] == 1)
	{{
		cI[idx] = cR[idx];
	}}
}}

__global__ void calcMedialPersistence(int *cR, int *cI, int *MPabs, float *MPrel)
{{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(idx < N && cR[idx] > 0)
	{{
		MPabs[idx] = cR[idx] - cI[idx];
		MPrel[idx] = 1. - cI[idx] / cR[idx];
	}}
}}

__global__ void markRemove(int *sFacets, int *sFaces, int *cI,
	int *MPabs, int *MPrel, int absT, float relT, int useMP, 
	int *rFacets, int *rFaces, int *sR)
{{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idx2;
	
	if(idx < N && sFacets[idx] == 1)
	{{
		if(useMP)
		{{
			int dim = dimension(idx);
			if(cI[idx] == 0 || MPabs[idx] < absT || MPrel[idx] < relT)
			{{
				rFacets[idx] = 1;
				sR[0] = 1;
			}}
		}}
		else
		{{
			rFacets[idx] = 1;
			sR[0] = 1;
		}}
		
		idx2 = idx + D1;
		if(idx2 < N && sFaces[idx2] == 1)
		{{
			if(useMP)
			{{
				int dim = dimension(idx2);
				if(cI[idx] == 0 || MPabs[idx] < absT || MPrel[idx] < relT)
				{{
					rFaces[idx2] = 1;
				}}
			}}
			else
			{{
				rFaces[idx2] = 1;
			}}
		}}
		idx2 = idx - D1;
		if(idx2 >= 0 && sFaces[idx2] == 1)
		{{
			if(useMP)
			{{
				int dim = dimension(idx2);
				if(cI[idx] == 0 || MPabs[idx] < absT || MPrel[idx] < relT)
				{{
					rFaces[idx2] = 1;
				}}
			}}
			else
			{{
				rFaces[idx2] = 1;
			}}
		}}
		
		idx2 = idx + D2;
		if(idx2 < N && sFaces[idx2] == 1)
		{{
			if(useMP)
			{{
				int dim = dimension(idx2);
				if(cI[idx] == 0 || MPabs[idx] < absT || MPrel[idx] < relT)
				{{
					rFaces[idx2] = 1;
				}}
			}}
			else
			{{
				rFaces[idx2] = 1;
			}}
		}}
		idx2 = idx - D2;
		if(idx2 >= 0 && sFaces[idx2] == 1)
		{{
			if(useMP)
			{{
				int dim = dimension(idx2);
				if(cI[idx] == 0 || MPabs[idx] < absT || MPrel[idx] < relT)
				{{
					rFaces[idx2] = 1;
				}}
			}}
			else
			{{
				rFaces[idx2] = 1;
			}}
		}}
		
		if(N3 > 1)
		{{
			idx2 = idx + D3;
			if(idx2 < N && sFaces[idx2] == 1)
			{{
				if(useMP)
				{{
					int dim = dimension(idx2);
					if(cI[idx] == 0 || MPabs[idx] < absT || MPrel[idx] < relT)
					{{
						rFaces[idx2] = 1;
					}}
				}}
				else
				{{
					rFaces[idx2] = 1;
				}}
			}}
			idx2 = idx - D3;
			if(idx2 >= 0 && sFaces[idx2] == 1)
			{{
				if(useMP)
				{{
					int dim = dimension(idx2);
					if(cI[idx] == 0 || MPabs[idx] < absT || MPrel[idx] < relT)
					{{
						rFaces[idx2] = 1;
					}}
				}}
				else
				{{
					rFaces[idx2] = 1;
				}}
			}}
		}}
	}}
}}

__global__ void doRemove(int *rFacets, int *rFaces, int *K)
{{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(idx < N && (rFacets[idx] > 0 || rFaces[idx] > 0))
	{{
		K[idx] = 0;
	}}
}}
