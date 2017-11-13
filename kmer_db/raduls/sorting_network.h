#pragma once
//#include "common.h"

#include "defs.h"
#include "record.h"
#include <vector>
#include <functional>
namespace raduls
{
	namespace small_sort
	{
		template<typename T, typename CompAndSwap>
		struct SortingNetwork
		{
			static const uint32 MAX_SUPPORTED_SIZE = 64;
			SortingNetwork();

			std::vector < std::function<void(T*)>> funcs;
			FORCE_INLINE void operator()(T* arr, uint32 size)
			{
				funcs[size](arr);
			}
			void sort2(T* arr);
			void sort3(T* arr);
			void sort4(T* arr);
			void sort5(T* arr);
			void sort6(T* arr);
			void sort7(T* arr);
			void sort8(T* arr);
			void sort9(T* arr);
			void sort10(T* arr);

			void sort11(T* arr);
			void sort12(T* arr);
			void sort13(T* arr);
			void sort14(T* arr);
			void sort15(T* arr);
			void sort16(T* arr);
			void sort17(T* arr);
			void sort18(T* arr);
			void sort19(T* arr);
			void sort20(T* arr);

			void sort21(T* arr);
			void sort22(T* arr);
			void sort23(T* arr);
			void sort24(T* arr);
			void sort25(T* arr);
			void sort26(T* arr);
			void sort27(T* arr);
			void sort28(T* arr);
			void sort29(T* arr);
			void sort30(T* arr);

			void sort31(T* arr);
			void sort32(T* arr);
			void sort33(T* arr);
			void sort34(T* arr);
			void sort35(T* arr);
			void sort36(T* arr);
			void sort37(T* arr);
			void sort38(T* arr);
			void sort39(T* arr);
			void sort40(T* arr);

			void sort41(T* arr);
			void sort42(T* arr);
			void sort43(T* arr);
			void sort44(T* arr);
			void sort45(T* arr);
			void sort46(T* arr);
			void sort47(T* arr);
			void sort48(T* arr);
			void sort49(T* arr);
			void sort50(T* arr);

			void sort51(T* arr);
			void sort52(T* arr);
			void sort53(T* arr);
			void sort54(T* arr);
			void sort55(T* arr);
			void sort56(T* arr);
			void sort57(T* arr);
			void sort58(T* arr);
			void sort59(T* arr);
			void sort60(T* arr);

			void sort61(T* arr);
			void sort62(T* arr);
			void sort63(T* arr);
			void sort64(T* arr);
		};
	}
}