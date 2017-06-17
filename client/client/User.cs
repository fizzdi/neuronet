namespace client
{
    public class User
    {
        private int _id;
        private string _login;
        private string _password;
        private string _teamname;

        public int ID { get { return _id; } }
        public string Login { get { return _login; } }
        public string Password { get { return _password; } }
        public string TeamName { get { return _teamname; } set { _teamname = value; } }

        public User(int id, string login, string password, string teamname)
        {
            _id = id;
            _login = login;
            _password = password;
            _teamname = teamname;
        }

        public bool EqualPassowrd(string password)
        {
            return password == _password;
        }
    }
}