namespace client
{
    partial class ShowSubmitCode
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.rtf_submitCode = new System.Windows.Forms.RichTextBox();
            this.b_selectAll = new System.Windows.Forms.Button();
            this.SuspendLayout();
            // 
            // rtf_submitCode
            // 
            this.rtf_submitCode.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.rtf_submitCode.BackColor = System.Drawing.SystemColors.ControlLightLight;
            this.rtf_submitCode.Location = new System.Drawing.Point(13, 13);
            this.rtf_submitCode.Name = "rtf_submitCode";
            this.rtf_submitCode.ReadOnly = true;
            this.rtf_submitCode.ShowSelectionMargin = true;
            this.rtf_submitCode.Size = new System.Drawing.Size(259, 207);
            this.rtf_submitCode.TabIndex = 0;
            this.rtf_submitCode.Text = "";
            // 
            // b_selectAll
            // 
            this.b_selectAll.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left)));
            this.b_selectAll.Location = new System.Drawing.Point(80, 226);
            this.b_selectAll.Name = "b_selectAll";
            this.b_selectAll.Size = new System.Drawing.Size(140, 23);
            this.b_selectAll.TabIndex = 1;
            this.b_selectAll.Text = "Скопировать в буффер";
            this.b_selectAll.UseVisualStyleBackColor = true;
            this.b_selectAll.Click += new System.EventHandler(this.b_selectAll_Click);
            // 
            // ShowSubmitCode
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(284, 262);
            this.Controls.Add(this.b_selectAll);
            this.Controls.Add(this.rtf_submitCode);
            this.MinimumSize = new System.Drawing.Size(300, 300);
            this.Name = "ShowSubmitCode";
            this.Load += new System.EventHandler(this.ShowSubmitCode_Load);
            this.Resize += new System.EventHandler(this.ShowSubmitCode_Resize);
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.RichTextBox rtf_submitCode;
        private System.Windows.Forms.Button b_selectAll;
    }
}