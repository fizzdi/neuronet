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

        World world = new World();
        List<List<object>> players = new List<List<object>>();
        List<bool> nums = new List<bool>(4);

        public MainForm()
        {
            InitializeComponent();
        }

        private void MainForm_Load(object sender, EventArgs e)
        {
            File.Delete("debug.txt");
            if (Directory.Exists("solutions"))
                Directory.Delete("solutions", true);
            if (Directory.Exists("players"))
                Directory.Delete("players", true);
            Directory.CreateDirectory("players");
            Directory.CreateDirectory("solutions");

            rb_solve_CheckedChanged(sender, e);
            for (int i = 0; i < 4; ++i)
                nums.Add(false);

        }

        private void b_localRun_Click(object sender, EventArgs e)
        {
            world.freePlayers();
            world.start();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            world.freePlayers();
            int n = 1;
            foreach (var pl in players)
            {
                if (pl.Count == 1)
                    File.Copy((string)pl[0], (string)pl[0]);
                else
                {
                    string pl_str = "player" + n;
                    Directory.CreateDirectory("solutions\\" + pl_str);
                    foreach (var file in Directory.GetFiles("solutions\\" + pl_str))
                        File.Delete(file);

                    foreach (var file in Directory.GetFiles("template"))
                        File.Copy(file, "solutions\\" + pl_str + "\\" + Path.GetFileName(file), true);
                    String res = File.ReadAllText("template\\MLP.h");
                    res = res.Replace("{0}", ((Decimal)pl[0]).ToString()); //SENSOR_COUNT
                    res = res.Replace("{1}", ((Decimal)pl[1]).ToString()); //HIDDEN_NEURON_COUNT
                    res = res.Replace("{2}", ((Decimal)pl[2]).ToString()); //HIDDEN_LAYER_COUNT
                    res = res.Replace("{3}", ((Decimal)pl[3]).ToString()); //TEST_COUNT
                    res = res.Replace("{4}", ((Decimal)pl[4]).ToString()); //MAX_TEST_COUNT
                    res = res.Replace("{5}", ((Decimal)pl[5]).ToString()); //TRAIN_EPOCH
                    res = res.Replace("{6}", ((Decimal)pl[6]).ToString()); //TRAIN_PERION
                    res = res.Replace("{7}", ((Decimal)pl[7]).ToString("0.00000000", CultureInfo.InvariantCulture)); //TRAIN_EPS
                    res = res.Replace("{8}", ((Decimal)pl[8]).ToString()); //END_TRAIN_TICK
                    res = res.Replace("{9}", ((int)pl[9]).ToString()); //FUN_ACT
                    res = res.Replace("{10}", ((Decimal)pl[10]).ToString("0.00000000", CultureInfo.InvariantCulture)); //RMS_GAMMA
                    res = res.Replace("{11}", ((Decimal)pl[11]).ToString("0.00000000", CultureInfo.InvariantCulture)); //RMS_LEARNRATE
                    res = res.Replace("{12}", ((Decimal)pl[12]).ToString("0.00000000", CultureInfo.InvariantCulture)); //RMS_EPSILON
                    File.WriteAllText("solutions\\" + pl_str + "\\MLP.h", res);

                    res = File.ReadAllText("template\\MyPlayer.cpp");
                    res = res.Replace("{0}", ((Decimal)pl[13]).ToString("0.00000000", CultureInfo.InvariantCulture));
                    res = res.Replace("{1}", (string)pl[14]);
                    res = res.Replace("{2}", (string)pl[15]);
                    File.WriteAllText("solutions\\" + pl_str + "\\MyPlayer.cpp", res);
                }
                n++;
            }

            //maybe debug
            world.freePlayers();
            foreach (var file in Directory.GetFiles("players", "*.dll"))
                File.Delete(file);
            //-----

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

        private void button2_Click(object sender, EventArgs e)
        {
            world.freePlayers();
            button1_Click(sender, e);
            world.start();
        }

        private void button3_Click(object sender, EventArgs e)
        {
            world.start();
        }

        private void rb_solve_CheckedChanged(object sender, EventArgs e)
        {

        }

        private void b_net_save_Click(object sender, EventArgs e)
        {



        }

        private void b_add_player_Click(object sender, EventArgs e)
        {
            Edit_player frm = new Edit_player();
            frm.ShowDialog();
            if (frm.is_net)
                players.Add(frm.p_opt);
            else
            {
                List<object> tmp = new List<object>();
                tmp.Add(frm.filename);
                players.Add(tmp);
            }
            for (int i = 0; i < 4; ++i)
            {
                if (nums[i]) continue;
                lb_players.Items.Add("player" + i);
                nums[i] = true;
                break;
            }
            if (lb_players.Items.Count == 4)
                b_add_player.Enabled = false;
        }

        private void b_edit_player_Click(object sender, EventArgs e)
        {
            Edit_player frm = new Edit_player(players[lb_players.SelectedIndex]);
            frm.ShowDialog();
            if (frm.is_net)
                players[lb_players.SelectedIndex] = frm.p_opt;
            else
            {
                List<object> tmp = new List<object>();
                tmp.Add(frm.filename);
                players[lb_players.SelectedIndex] = tmp;
            }
        }

        private void b_remove_player_Click(object sender, EventArgs e)
        {
            string str = lb_players.SelectedItem.ToString();
            int num = str[str.Length - 1] - '0';
            nums[num] = false;
            players.RemoveAt(lb_players.SelectedIndex);
            b_add_player.Enabled = true;
            lb_players.Items.RemoveAt(lb_players.SelectedIndex);
        }

        private void lb_players_SelectedIndexChanged(object sender, EventArgs e)
        {
            b_edit_player.Enabled = b_remove_player.Enabled = lb_players.SelectedIndex != -1;
        }
    }
}
