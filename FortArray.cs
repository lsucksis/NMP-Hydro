/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2021/5/5
 * Time: 14:41
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;

namespace NoahMP
{
	/// <summary>
	/// 模仿fortran版本的数组，索引从1起始，或者从一负数起始
	/// Description of FortArray.
	/// </summary>
	public class FortDoubleArray
	{
		public double[] data = null;
		int start = 1;
		public FortDoubleArray(int start, int end)
		{
			data = new double[end - start + 1];
			this.start = start;
		}
		public FortDoubleArray(double[] values)
		{
			data = values;
			this.start = 1;
		}
		public double this[int index] {
			get {
				return data[index - start];
			}
			set {
				data[index - start] = value;
			}			
		}
		public int Count {
			get {
				return data.Length;
			}
		}
	}
	public class FortIntArray
	{
		public int[] data = null;
		int start = 1;
		public FortIntArray(int start, int end)
		{
			data = new int[end - start + 1];
			this.start = start;
		}
		public int this[int index] {
			get {
				return data[index - start];
			}
			set {
				data[index - start] = value;
			}			
		}
		public int Count {
			get {
				return data.Length;
			}
		}
	}
}
