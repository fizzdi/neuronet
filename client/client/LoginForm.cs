using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace client
{
    public partial class LoginForm : Form
    {
        private Server server;
        private User user = null;
        private bool isEnter = false;

        public LoginForm(Server server)
        {
            this.server = server;
            InitializeComponent();
        }

        private void b_ok_Click(object sender, EventArgs e)
        {
            if (tb_login.Text == "")
            {
                MessageBox.Show("Введите логин", "Ошибка входа", MessageBoxButtons.OK, MessageBoxIcon.Information);
                tb_login.Focus();
                return;
            }

            if (tb_password.Text == "")
            {
                MessageBox.Show("Введите пароль", "Ошибка входа", MessageBoxButtons.OK, MessageBoxIcon.Information);
                tb_password.Focus();
                return;
            }

            if (!verifLogin(tb_login.Text))
            {
                MessageBox.Show("Неверный логин!", "Ошибка входа", MessageBoxButtons.OK, MessageBoxIcon.Error);
                tb_login.Focus();
                return;
            }

            
            try
            {
                 user = server.getUser(tb_login.Text);
            }
            catch(Exception ex)
            {
                Logger.LogMessage(ex);
                MessageBox.Show("Неверный логин!", "Ошибка входа", MessageBoxButtons.OK, MessageBoxIcon.Error);
                tb_login.Focus();
                tb_password.Clear();
                return;
            }

            if (!verifPassword(tb_password.Text) || !user.EqualPassowrd(tb_password.Text))
            {
                MessageBox.Show("Неверный пароль", "Ошибка входа", MessageBoxButtons.OK, MessageBoxIcon.Error);
                tb_password.Clear();
                tb_password.Focus();
                return;
            }

            isEnter = true;
            (this.Owner as MainForm).currentUser = user;
            this.Close();
        }

        private bool verifPassword(string password)
        {
            for (int i = 0; i < password.Length; ++i)
                if (!Char.IsLetterOrDigit(password[i]))
                    return false;
            return true;
        }

        private bool verifLogin(string login)
        {
            for (int i = 0; i < login.Length; ++i)
                if (!Char.IsLetterOrDigit(login[i]))
                    return false;
            return true;
        }

        private void LoginForm_FormClosed(object sender, FormClosedEventArgs e)
        {
            if (!isEnter)
            {
                Owner.Close();
                this.Close();
            }
        }
    }
}
