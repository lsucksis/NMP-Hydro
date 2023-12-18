/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2016/11/26
 * Time: 20:57
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;

namespace NoahMP
{
	class MyDriver:Driver
	{
		public override void InitModelCells2()
		{
			Cells = new GridCellReal[NRows, NCols];
			for (int row = 0; row < NRows; row++) {
				for (int col = 0; col < NCols; col++) {
					Cells[row, col] = new GridCellReal();
				}
			}
		}
	}
	class Program
	{
		public static void Main(string[] args)
		{
			Console.WriteLine("=================================================================================================");
			Console.WriteLine("Welcome! This model is translated from the original fortran code of Noah-MP by Yonghe Liu, at Henan Polytechnic University.");
			Console.WriteLine("欢迎您使用。本模型是由河南理工大学刘永和在WRF-Hydro 3.0基础上翻译而来。\n" +
			                  "---2018年1-2月,由Noah-MP的fortran代码基础上直接翻译而成，并经过了严格的调试。当时个人所承担的教学任务已结束，为了弄清楚代码逻辑，不分昼夜共耗时1月有余。\n" +
			"---2019年时调整了数组的结构(使原来Fortran中以任何下标起始的方式与C#中的数组下标变得一致)，造成几乎重新翻译了一遍，但是这次调整后，代码开始更适合调试，因而调试工作量更小了，工作进度加快较多。最终发现，原来2018年的一个多月" +
			"的工作基本上是失败的。\n" +
			"---2021年时又进行了大量调试，以保证此Noah-MP与原Fortran版之间运行结果的一致性。\n" +
			"---2022年4月以及8-11月又做了调试，又发现了不少代码错误，其间还进行了多次重要的调整；解决了重启续模拟不衔接问题，即重启续模拟与一直连续运行具有完全" +
			"相同的效果。\n" +
			"---2022年11月18日，发现了一年以来一直寻找的bug，终于使C#版与Fortran版的运行结果几乎一致了。与Noah-MP原版运行结果不同，很大程度上是由于数值误差以及影响较微小的bug与造成。\n" +
			"---2022年12月，开发了CarbonSink模块，用于模拟岩溶碳汇。同时修改了动态植被模块，更新为WRF-Hydro5.2的植被模块版本。\n" +
			"---2022年6月22日（夏至日下午）开发了打包版。");
			Console.WriteLine("=================================================================================================");
			MyDriver drv = new MyDriver();			
			drv.InitStatus();
			drv.InitModelCells2();
			drv.RunModel();
			
			// TODO: Implement Functionality Here
			
			Console.Write("The simulation completed . . . ");
			
			
			//下面是河网汇流的代码
			//RiverRouting model=new RiverRouting();
			
		}
	}
}