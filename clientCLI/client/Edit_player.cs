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

namespace client
{
    public partial class Edit_player : Form
    {
        public List<object> p_opt = new List<object>();
        public Edit_player()
        {
            InitializeComponent();
            cb_train_type.SelectedIndex = 0;
            cb_fun_act.SelectedIndex = 0;
        }

        public Edit_player(List<object> p_opt)
        {
            InitializeComponent();
            if (p_opt.Count <= 2)
            {
                string filename = (string)p_opt[0];
                rb_solve.Checked = true;
                try
                {
                    l_filename.Text = (new FileInfo(filename)).Name;
                    using (StreamReader reader = new StreamReader(File.Open(filename, FileMode.Open)))
                        rtf_solve.Text = reader.ReadToEnd();
                }
                catch (Exception ex)
                {
                    Logger.LogMessage(ex);
                }
            }
            else
            {
                this.p_opt = p_opt;
                nud_sensors.Value = (Decimal)p_opt[0];
                nud_hid_neurons.Value = (Decimal)p_opt[1];
                nud_hid_layers.Value = (Decimal)p_opt[2];
                nud_test_count.Value = (Decimal)p_opt[3];
                nud_max_test_count.Value = (Decimal)p_opt[4];
                nud_train_epoch.Value = (Decimal)p_opt[5];
                nud_train_period.Value = (Decimal)p_opt[6];
                nud_train_eps.Value = (Decimal)p_opt[7];
                nud_end_train.Value = (Decimal)p_opt[8];
                cb_fun_act.SelectedIndex = p_opt[9].GetType() == typeof(Decimal) ? Decimal.ToInt32((Decimal)p_opt[9]) : (int)p_opt[9];
                nud_rms_gamma.Value = (Decimal)p_opt[10];
                nud_rms_learnrate.Value = (Decimal)p_opt[11];
                nud_rms_eps.Value = (Decimal)p_opt[12];
                nud_q_learn.Value = (Decimal)p_opt[13];
                tb_name.Text = (String)p_opt[14];
                cb_train_type.SelectedIndex = cb_train_type.Items.IndexOf(p_opt[15]);
                rb_net.Checked = true;
            }
        }

        private void radioButton1_CheckedChanged(object sender, EventArgs e)
        {
            nud_rms_gamma.Enabled = cb_train_type.SelectedIndex == 0;
            nud_rms_learnrate.Enabled = cb_train_type.SelectedIndex == 0;
            nud_rms_eps.Enabled = cb_train_type.SelectedIndex == 0;
            if (rb_solve.Checked)
            {
                p_solve.Visible = true;
                p_net.Visible = false;
            }
            else
            {
                p_net.Visible = true;
                p_solve.Visible = false;
            }
        }

        private void b_clear_Click(object sender, EventArgs e)
        {
            rtf_solve.Clear();
            l_filename.Text = "Файл не выбран";
        }

        private void b_open_Click(object sender, EventArgs e)
        {
            ofd_sourseFile.ShowDialog();
            try
            {
                l_filename.Text = (new FileInfo(ofd_sourseFile.FileName)).Name;
                using (StreamReader reader = new StreamReader(ofd_sourseFile.OpenFile()))
                    rtf_solve.Text = reader.ReadToEnd();
            }
            catch (Exception ex)
            {
                Logger.LogMessage(ex);
            }
        }

        private void b_save_Click(object sender, EventArgs e)
        {
            if (rb_net.Checked)
            {
                while (p_opt.Count < 16)
                    p_opt.Add(0);

                p_opt[0] = nud_sensors.Value;
                p_opt[1] = nud_hid_neurons.Value;
                p_opt[2] = nud_hid_layers.Value;
                p_opt[3] = nud_test_count.Value;
                p_opt[4] = nud_max_test_count.Value;
                p_opt[5] = nud_train_epoch.Value;
                p_opt[6] = nud_train_period.Value;
                p_opt[7] = nud_train_eps.Value;
                p_opt[8] = nud_end_train.Value;
                p_opt[9] = cb_fun_act.SelectedIndex;
                p_opt[10] = nud_rms_gamma.Value;
                p_opt[11] = nud_rms_learnrate.Value;
                p_opt[12] = nud_rms_eps.Value;
                p_opt[13] = nud_q_learn.Value;
                p_opt[14] = tb_name.Text;
                p_opt[15] = cb_train_type.SelectedItem;
            }
            else
            {

                p_opt.Clear();
                string filename;
                Directory.CreateDirectory("temp");
                if (File.Exists("temp\\" + l_filename.Text))
                    filename = "temp\\" + l_filename.Text;
                else
                    filename = string.Format("temp\\{0}{1}.cpp", DateTime.Now.Millisecond.ToString(), DateTime.Now.Minute.ToString());
                File.WriteAllText(filename, rtf_solve.Text);
                p_opt.Add(filename);
            }
            DialogResult = DialogResult.OK;
            Close();
        }

        private void Edit_player_Load(object sender, EventArgs e)
        {
        }

        private void cb_train_type_SelectedIndexChanged(object sender, EventArgs e)
        {
            nud_rms_eps.Enabled = nud_rms_gamma.Enabled = nud_rms_learnrate.Enabled = cb_train_type.SelectedIndex == 0;
        }
    }
}
