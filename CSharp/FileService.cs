using Newtonsoft.Json;
using System.IO;

namespace CSharp
{
    static class FileService
    {
        public static void SaveInJson(string path, object value)
        {
            string json = JsonConvert.SerializeObject(value);
            File.WriteAllText(path, json);
        }

        public static T ReadFromJson<T>(string path)
        {
            string json = File.ReadAllText(path);
            T value = JsonConvert.DeserializeObject<T>(json);
            return value;
        }
    }
}
