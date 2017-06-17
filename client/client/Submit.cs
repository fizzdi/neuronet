namespace client
{
    public class Submit
    {
        private int _number;
        private string _status;
        private string _code;

        public int Number { get { return _number; } }
        public string Status { get { return _status; } }
        public string Code { get { return _code; } }


        public Submit(int number, string status, string code)
        {
            _number = number;
            _status = status;
            _code = code;
        }
    }
}