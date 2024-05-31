#pragma once
#include <cstring>
#include <cstdint>

#ifdef _MSC_VER  
#include <intrin.h>
#endif 

namespace InstrSetDetect
{
	enum class Instr { NotSet, SSE, SSE2, SSE3, SSE4_1, SSE4_2, AVX, AVX2, NEON };

#ifdef ARCH_X64
	static void cpuid(uint32_t eax, uint32_t ecx, uint32_t* abcd)
	{
#if _MSC_VER
		__cpuidex((int*)abcd, eax, ecx);
#else
		uint32_t ebx=0, edx;
		__asm__("cpuid" : "+b" (ebx), "+a" (eax), "+c" (ecx), "=d" (edx));
		abcd[0] = eax; abcd[1] = ebx; abcd[2] = ecx; abcd[3] = edx;
#endif
	}

	static Instr GetInstr()
	{
		Instr instr{ Instr::NotSet };
		uint32_t abcd[4]{};
		cpuid(1, 0, abcd);

		uint32_t edx = abcd[3];

		uint32_t ecx = abcd[2];

		if (((edx >> 25) & 1) == 0) return instr;
		instr = Instr::SSE;

		if (((edx >> 26) & 1) == 0) return instr;
		instr = Instr::SSE2;

		if (((ecx >> 0) & 1) == 0) return instr;
		instr = Instr::SSE3;

		if (((ecx >> 19) & 1) == 0) return instr;
		instr = Instr::SSE4_1;

		if (((ecx >> 20) & 1) == 0) return instr;
		instr = Instr::SSE4_2;

		if (((ecx >> 28) & 1) == 0) return instr;
		instr = Instr::AVX;

		cpuid(7, 0, abcd);

		uint32_t ebx = abcd[1];
		if (((ebx >> 5) & 1) == 0) return instr;
		instr = Instr::AVX2;

		return instr;
	}
#else
	static Instr GetInstr()
	{
		return Instr::NEON;
	}
#endif
};
