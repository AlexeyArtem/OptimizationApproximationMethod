using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Text.Json;

namespace CSharp
{
    class Program
    {
        static void Main(string[] args)
        {
            List item = await JsonFileReader.ReadAsync<Item>(@"C:\myFile.json");

        }
    }

    public static class JsonFileReader
    {
        public static async Task<T> ReadAsync<T>(string filePath)
        {
            FileStream stream = File.OpenRead(filePath);
            return await JsonSerializer.DeserializeAsync<T>(stream);
        }
    }
}
