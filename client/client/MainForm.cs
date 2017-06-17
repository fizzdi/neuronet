using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Globalization;
using System.Runtime.InteropServices;

namespace client
{
    public partial class MainForm : Form
    {
        //[DllImport("kernel32.dll", EntryPoint = "LoadLibrary", SetLastError = true)]
        [DllImport("kernel32.dll")]
        public static extern IntPtr LoadLibrary(string dllToLoad);
        [DllImport("kernel32.dll")]
        public static extern bool FreeLibrary(IntPtr hModule);

        Server server;
        World world = new World();
        bool[] isBotStrategy = new bool[4];
        private string ip;
        ServerConfig cfg;
        public User currentUser;
        Dictionary<string, int> Bots = new Dictionary<string, int>();
        Dictionary<int, Submit> submits;

        public MainForm()
        {
            InitializeComponent();
        }

        private void MainForm_Load(object sender, EventArgs e)
        {
            try
            {
                cfg = new ServerConfig();
            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.Message, "Ошибка инициализации", MessageBoxButtons.OK, MessageBoxIcon.Error);
                this.Close();
                Application.Exit();
                return;
            }

            server = new Server(cfg);
            //server.connect();
            // server.addBots(lb_local_bots);
            world.server = server;
        }

        private void b_downloadFile_Click(object sender, EventArgs e)
        {
            ofd_sourseFile.ShowDialog();
            try
            {
                l_sourceFile.Text = (new FileInfo(ofd_sourseFile.FileName)).Name;
                using (StreamReader reader = new StreamReader(ofd_sourseFile.OpenFile()))
                    rtf_sourse.Text = reader.ReadToEnd();
            }
            catch (Exception ex)
            {
                Logger.LogMessage(ex);
            }
        }

        private void b_send_Click(object sender, EventArgs e)
        {
            if (rtf_sourse.Text == "")
            {
                MessageBox.Show("Введите код решения", "Ошибка отправления решения", MessageBoxButtons.OK, MessageBoxIcon.Exclamation);
                return;
            }

            server.uploadSubmit(1, rtf_sourse.Text);
            MessageBox.Show("Решение отправлено успешно", "Отправка решения", MessageBoxButtons.OK, MessageBoxIcon.Information);
        }

        private void b_clear_send_Click(object sender, EventArgs e)
        {
            rtf_sourse.Clear();
            l_sourceFile.Text = "Файл не выбран";
        }

        private void b_reload_history_Click(object sender, EventArgs e)
        {
            dgv_history.Rows.Clear();
            submits = server.getSubmits(currentUser.ID);
            foreach (var cur in submits)
                dgv_history.Rows.Add(dgv_history.RowCount + 1, cur.Value.Number, cur.Value.Status, "Посмотреть");
        }

        private void dgv_history_Resize(object sender, EventArgs e)
        {
            col_submitState.Width = dgv_history.Width - col_numberOfSubmit.Width - col_timeOfSubmit.Width - col_showButton.Width - 3;
        }

        private void dgv_history_CellContentClick(object sender, DataGridViewCellEventArgs e)
        {
            if (e.RowIndex == -1)
                return;

            if (e.ColumnIndex == col_showButton.Index)
            {
                ShowSubmitCode ssc = new ShowSubmitCode(submits[(int)dgv_history.Rows[e.RowIndex].Cells[col_numberOfSubmit.Index].Value]);
                ssc.ShowDialog(this);
            }
        }

        private void tab_local_run_Resize(object sender, EventArgs e)
        {
            gb_local_myStrategy.Height = tab_local_run.Height / 2;
            gb_local_bots.Height = tab_local_run.Height / 2 - gb_local_myStrategy.Top;
            gb_local_bots.Top = gb_local_myStrategy.Height + gb_local_myStrategy.Top;
            int oneBlock = tab_local_run.Width / 7;
            gb_local_myStrategy.Width = gb_local_bots.Width = oneBlock * 3 - gb_local_bots.Left;
            b_local_left.Left = b_local_right.Left = b_localRun.Left = gb_local_bots.Left + gb_local_bots.Width + 3;
            b_local_left.Width = b_local_right.Width = b_localRun.Width = oneBlock - 6;
            gb_local_selectedStrategy.Left = b_local_left.Left + b_local_left.Width + 3;
            gb_local_selectedStrategy.Width = tab_local_run.Width - gb_local_bots.Left - gb_local_selectedStrategy.Left;
        }

        private void b_localRight_Click(object sender, EventArgs e)
        {
          /*  if (lb_local_bots.SelectedIndex == -1)
            {
                //Добавление стратегии игрока
                lb_local_selectedStrategy.Items.Add(lb_local_myStrategy.SelectedItem);
                int curSelect = lb_local_myStrategy.SelectedIndex;
                lb_local_myStrategy.SelectedIndex = Math.Min(curSelect, lb_local_myStrategy.Items.Count - 1);
            }
            else
            {
                //Добавление стратегии ботов
                lb_local_selectedStrategy.Items.Add(lb_local_bots.SelectedItem);
                int curSelect = lb_local_bots.SelectedIndex;
                lb_local_bots.SelectedIndex = Math.Min(curSelect, lb_local_bots.Items.Count - 1);
            }

            if (lb_local_selectedStrategy.Items.Count > 0)
            {
                b_localRun.Enabled = true;
            }

            if (lb_local_selectedStrategy.Items.Count == 4)
            {
                b_local_right.Enabled = false;
            }*/
        }

        private void lb_myStrategy_Click(object sender, EventArgs e)
        {
            lb_local_bots.SelectedIndex = -1;
        }

        private void lb_bots_Click(object sender, EventArgs e)
        {
            lb_local_myStrategy.SelectedIndex = -1;
        }

        private void b_localLeft_Click(object sender, EventArgs e)
        {
            /*int curSelected = lb_local_selectedStrategy.SelectedIndex;

            lb_local_selectedStrategy.Items.RemoveAt(lb_local_selectedStrategy.SelectedIndex);

            if (lb_local_selectedStrategy.Items.Count == 0)
            {
                b_localRun.Enabled = false;
                b_local_left.Enabled = false;
            }
            else
                lb_local_selectedStrategy.SelectedIndex = Math.Min(curSelected, lb_local_selectedStrategy.Items.Count - 1);*/
        }

        private void lb_myStrategy_SelectedIndexChanged(object sender, EventArgs e)
        {
            //b_local_right.Enabled = lb_local_selectedStrategy.Items.Count != 4 && !(lb_local_myStrategy.SelectedIndex == -1 && lb_local_bots.SelectedIndex == lb_local_myStrategy.SelectedIndex);
        }

        private void lb_bots_SelectedIndexChanged(object sender, EventArgs e)
        {
            //b_local_right.Enabled = lb_local_selectedStrategy.Items.Count != 4 && !(lb_local_bots.SelectedIndex == -1 && lb_local_bots.SelectedIndex == lb_local_myStrategy.SelectedIndex);
        }

        private void lb_selectedStrategy_SelectedIndexChanged(object sender, EventArgs e)
        {
           // b_local_left.Enabled = lb_local_selectedStrategy.SelectedIndex != -1;
          //  lb_local_myStrategy.SelectedIndex = lb_local_bots.SelectedIndex = -1;
        }

        private void b_localRun_Click(object sender, EventArgs e)
        {
            /*foreach (string cur in lb_local_selectedStrategy.Items)
            {
                if (!Char.IsDigit(cur[0]))
                    world.addBots(cur);
                else
                {
                    world.addPlayers(int.Parse(cur));
                }
            }*/
            world.start();
        }

        private void tab_global_run_Resize(object sender, EventArgs e)
        {
            gb_global_myStrategy.Height = tab_global_run.Height / 2;
            gb_global_otherStrategy.Height = tab_global_run.Height / 2 - gb_global_myStrategy.Top;
            gb_global_otherStrategy.Top = gb_global_myStrategy.Height + gb_global_myStrategy.Top;
            int oneBlock = tab_global_run.Width / 7;
            gb_global_myStrategy.Width = gb_global_otherStrategy.Width = oneBlock * 3 - gb_global_otherStrategy.Left;
            b_global_left.Left = b_global_right.Left = b_globalRun.Left = gb_global_otherStrategy.Left + gb_global_otherStrategy.Width + 3;
            b_global_left.Width = b_global_right.Width = b_globalRun.Width = oneBlock - 6;
            gb_global_selectedStrategy.Left = b_global_left.Left + b_global_left.Width + 3;
            gb_global_selectedStrategy.Width = tab_global_run.Width - gb_global_otherStrategy.Left - gb_global_selectedStrategy.Left;
        }

        private void lb_global_myStrategy_Click(object sender, EventArgs e)
        {
            lb_global_otherStrategy.SelectedIndex = -1;
        }

        private void lb_global_otherStrategy_Click(object sender, EventArgs e)
        {
            lb_global_myStrategy.SelectedIndex = -1;
        }

        private void b_global_right_Click(object sender, EventArgs e)
        {
            if (lb_global_otherStrategy.SelectedIndex == -1)
            {
                //Добавление стратегии игрока
                lb_global_selectedStrategy.Items.Add(lb_global_myStrategy.SelectedItem);
                int curSelect = lb_global_myStrategy.SelectedIndex;
                lb_global_myStrategy.SelectedIndex = Math.Min(curSelect, lb_global_myStrategy.Items.Count - 1);
            }
            else
            {
                //Добавление стратегии других участников
                lb_global_selectedStrategy.Items.Add(lb_global_otherStrategy.SelectedItem);
                int curSelect = lb_global_otherStrategy.SelectedIndex;
                lb_global_otherStrategy.SelectedIndex = Math.Min(curSelect, lb_global_otherStrategy.Items.Count - 1);
            }

            if (lb_global_selectedStrategy.Items.Count > 0)
            {
                b_globalRun.Enabled = true;
            }

            if (lb_global_selectedStrategy.Items.Count == 4)
            {
                //b_global_right.Enabled = false;
            }
        }

        private void b_global_left_Click(object sender, EventArgs e)
        {
            int curSelected = lb_global_selectedStrategy.SelectedIndex;

            lb_global_selectedStrategy.Items.RemoveAt(lb_global_selectedStrategy.SelectedIndex);

            if (lb_global_selectedStrategy.Items.Count == 0)
            {
                b_globalRun.Enabled = false;
                b_global_left.Enabled = false;
            }
            else
                lb_global_selectedStrategy.SelectedIndex = Math.Min(curSelected, lb_global_selectedStrategy.Items.Count - 1);
        }

        private void lb_global_myStrategy_SelectedIndexChanged(object sender, EventArgs e)
        {
            b_global_right.Enabled = lb_global_selectedStrategy.Items.Count != 4 && !(lb_global_myStrategy.SelectedIndex == -1 && lb_global_otherStrategy.SelectedIndex == lb_global_myStrategy.SelectedIndex);
        }

        private void lb_global_otherStrategy_SelectedIndexChanged(object sender, EventArgs e)
        {
            b_global_right.Enabled = lb_global_selectedStrategy.Items.Count != 4 && !(lb_global_otherStrategy.SelectedIndex == -1 && lb_global_otherStrategy.SelectedIndex == lb_global_myStrategy.SelectedIndex);
        }

        private void lb_global_selectedStrategy_SelectedIndexChanged(object sender, EventArgs e)
        {
            b_global_left.Enabled = lb_global_selectedStrategy.SelectedIndex != -1;
            lb_global_myStrategy.SelectedIndex = lb_global_otherStrategy.SelectedIndex = -1;
        }

        private void b_globalRun_Click(object sender, EventArgs e)
        {
            List<int> submits = new List<int>();
            foreach (string cur in lb_global_selectedStrategy.Items)
            {
                world.addPlayers(int.Parse(cur.Substring(cur.IndexOf('#') + 1)));
            }
            world.start();
        }

        private void tabcontrol_SelectedIndexChanged(object sender, EventArgs e)
        {
           /* if (tabcontrol.SelectedTab == tab_local_run)
            {
                lb_local_bots.Items.Clear();
                server.addBots(lb_local_bots);
                lb_local_myStrategy.Items.Clear();
                lb_local_myStrategy.Items.AddRange(server.getSubmits(currentUser.ID).Values.Where(t => t.Status != "Compile Error").OrderByDescending(t => t.Number).Select(t => t.Number.ToString()).ToArray());
            }

            if (tabcontrol.SelectedTab == tab_global_run)
            {
                lb_global_otherStrategy.Items.Clear();
                server.addOtherPlayersSubmits(currentUser.ID, lb_global_otherStrategy);
                lb_global_myStrategy.Items.Clear();
                lb_global_myStrategy.Items.AddRange(server.getSubmits(currentUser.ID).Values.Where(t => t.Status != "Compile Error").OrderByDescending(t => t.Number).Select(t => currentUser.TeamName + "#" + t.Number.ToString()).ToArray());
            }

            if (tabcontrol.SelectedTab == tab_history)
            {
                b_reload_history.PerformClick();
            }

            if (tabcontrol.SelectedTab == tab_settings)
            {
                tb_GroupName.Text = currentUser.TeamName;
            }*/
        }

        private void b_save_settings_Click(object sender, EventArgs e)
        {
            if (tb_GroupName.Text == "")
            {
                MessageBox.Show("Введите название команды", "Ошибка сохранения настроек", MessageBoxButtons.OK, MessageBoxIcon.Warning);
                return;
            }

            if (tb_GroupName.Text.IndexOf('#') != -1)
            {
                MessageBox.Show("Недопустимый символ в названии команды", "Ошибка сохранения настроек", MessageBoxButtons.OK, MessageBoxIcon.Error);
                return;
            }

            server.ChangeTeamName(currentUser, tb_GroupName.Text);
        }


        private void button1_Click(object sender, EventArgs e)
        {
            var lib = LoadLibrary("solve_compiller.dll");

            using (StreamWriter wr = new StreamWriter("dlldebug.txt"))
            {
                if (lib == IntPtr.Zero)
                    wr.WriteLine(Marshal.GetLastWin32Error()); // Error Code while loading DLL
                else
                    wr.WriteLine("Yes " + lib);  // Loading done !
            }
                if (lib != IntPtr.Zero)
                    FreeLibrary(lib);
        }
    }
}
