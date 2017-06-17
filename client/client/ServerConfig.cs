using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net;
using System.Text;
using System.Threading.Tasks;

namespace client
{
    public class ServerConfig
    {
        private IPAddress _ip;
        private string _username;
        private string _password;
        private string _bd_name;

        public IPAddress IP { get { return _ip; } }
        public string username { get { return _username; } }
        public string password { get { return _password; } }
        public string bdName { get { return _bd_name; } }

        public ServerConfig(string configFile = "client.cfg")
        {
            using (StreamReader reader = new StreamReader(configFile))
            {
                _ip = IPAddress.Parse(reader.ReadLine());
                _username = reader.ReadLine();
                _password = reader.ReadLine();
                _bd_name = reader.ReadLine();
            }
            if (_ip == null || _username == null || _password == null || _bd_name == null)
                throw new Exception("Неверные настройки клиента");
        }
    }
}
