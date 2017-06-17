using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;
using System.Text.RegularExpressions;

namespace client
{
    public class Server
    {
        ServerConfig cfg;

        public Server(ServerConfig cfg)
        {
            this.cfg = cfg;
        }

        public bool connect()
        {
            try
            {
               /* MySqlConnectionStringBuilder mcsb = new MySqlConnectionStringBuilder();
                mcsb.Server = cfg.IP.ToString();
                mcsb.UserID = cfg.username;
                mcsb.Password = cfg.password;
                mcsb.Database = cfg.bdName;
                mcsb.CharacterSet = "utf8";
                string connectionString = mcsb.ToString();
                conn = new MySqlConnection(mcsb.ToString());
                conn.Open();
                cmd = new MySqlCommand("", conn);*/
            }
            catch (Exception ex)
            {
                Logger.LogMessage(ex);
                MessageBox.Show(ex.Message, "Критическая ошибка", MessageBoxButtons.OK, MessageBoxIcon.Error);

                Application.ExitThread();
                return false;
            }
            return true;
        }

        public User getUser(string login)
        {
            /*if (reader != null && !reader.IsClosed)
            {
                reader.Close();
                reader = null;
            }
            cmd.CommandText = "select * from users where name='" + login + "'";
            reader = cmd.ExecuteReader();
           // reader.Read();
            User user = new User(reader.GetInt32("id"), reader.GetString("Name"), reader.GetString("Password"), reader.GetString("teamname"));
            reader.Close();
            return user;*/
            return new User(1,"","","");
        }

        public void downloadDll(int numberOfSubmit)
        {
           /* if (reader != null && !reader.IsClosed)
                reader.Close();
            cmd.CommandText = "select dll from submits where id='" + numberOfSubmit + "'";
            reader = cmd.ExecuteReader();

            if (!Directory.Exists("players"))
                Directory.CreateDirectory("players");

            FileStream dll = File.Create("players\\p" + Directory.GetFiles("players").Count() + ".dll");

            reader.Read();
            byte[] buff = (byte[])reader.GetValue(0);


            using (BinaryWriter writer = new BinaryWriter(dll))
            {
                writer.Write(buff);
            }
            reader.Close();*/
        }

        public void downloadBotDll(string name)
        {
            /*if (reader != null && !reader.IsClosed)
                reader.Close();
            cmd.CommandText = "select dll from bots where name='" + name + "'";
            reader = cmd.ExecuteReader();

            if (!Directory.Exists("players"))
                Directory.CreateDirectory("players");

            FileStream dll = File.Create("players\\b" + Directory.GetFiles("players").Count() + ".dll");

            reader.Read();
            byte[] buff = (byte[])reader.GetValue(0);

            using (BinaryWriter writer = new BinaryWriter(dll))
            {
                writer.Write(buff);
            }
            reader.Close();*/
        }

        public void uploadSubmit(int idUser, string code)
        {
           /* cmd.CommandText = "insert into submits (iduser, code) values('" + idUser + "','" + Regex.Escape(code).Replace("\'", "\\'") + "')";
            cmd.ExecuteNonQuery();*/
        }

        public Dictionary<int,Submit> getSubmits(int userId)
        {
            /*Dictionary<int, Submit> submits = new Dictionary<int, Submit>();
            
            if (reader != null && !reader.IsClosed)
                reader.Close();
            cmd.CommandText = "select submits.id, states.state, submits.code from submits, states where states.id = submits.compile and iduser='" + userId + "'";
            reader = cmd.ExecuteReader();

            while (reader.Read())
                submits.Add(reader.GetInt32("id"), new Submit(reader.GetInt32("id"), reader.GetString("state"), reader.GetString("code")));

            reader.Close();
            return submits;*/
            return new Dictionary<int, Submit>();
        } 

        public void addBots(ListBox lb)
        {
           /* if (reader != null && !reader.IsClosed)
                reader.Close();
            cmd.CommandText = "select name from bots order by name asc";
            reader = cmd.ExecuteReader();

            while (reader.Read())
                lb.Items.Add(reader.GetString(0));
            reader.Close();*/
        }

        public void addOtherPlayersSubmits(int iduser, ListBox lb)
        {
            /*if (reader != null && !reader.IsClosed)
                reader.Close();
            cmd.CommandText = "select concat(users.teamname, '#',s.id) submit from (select iduser, id from submits where compile > 0 and iduser <> '"+iduser+"'order by id desc) s, users where s.iduser = users.id group by s.iduser";
            reader = cmd.ExecuteReader();

            while (reader.Read())
                lb.Items.Add(reader.GetString(0));
            reader.Close();*/
        }

        public void ChangeTeamName(User user, string newTeamName)
        {
           /* cmd.CommandText = "update users set teamname='"+newTeamName+"' where id='"+user.ID+"'";
            cmd.ExecuteNonQuery();
            user.TeamName = newTeamName;*/
        }
    }
}
