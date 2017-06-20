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
    public partial class ShowSubmitCode : Form
    {
        private Submit submit;
        public ShowSubmitCode(Submit submit)
        {
            InitializeComponent();
            this.submit = submit;
        }

        private void ShowSubmitCode_Load(object sender, EventArgs e)
        {
            this.Icon = Owner.Icon;
            b_selectAll.Left = Width / 2 - b_selectAll.Width / 2;
            rtf_submitCode.Text = submit.Code;
            this.Text = "Просмотр решения №" + submit.Number + "(" + submit.Status + ")";
        }

        private void ShowSubmitCode_Resize(object sender, EventArgs e)
        {
            b_selectAll.Left = Width / 2 - b_selectAll.Width / 2;
        }

        private void b_selectAll_Click(object sender, EventArgs e)
        {
            rtf_submitCode.SelectAll();
            rtf_submitCode.Copy();
        }
    }
}
