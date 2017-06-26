namespace client
{
    partial class MainForm
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
            System.ComponentModel.ComponentResourceManager resources = new System.ComponentModel.ComponentResourceManager(typeof(MainForm));
            this.button2 = new System.Windows.Forms.Button();
            this.groupBox6 = new System.Windows.Forms.GroupBox();
            this.b_remove_player = new System.Windows.Forms.Button();
            this.b_edit_player = new System.Windows.Forms.Button();
            this.b_add_player = new System.Windows.Forms.Button();
            this.lb_players = new System.Windows.Forms.ListBox();
            this.groupBox6.SuspendLayout();
            this.SuspendLayout();
            // 
            // button2
            // 
            this.button2.Location = new System.Drawing.Point(189, 13);
            this.button2.Margin = new System.Windows.Forms.Padding(4);
            this.button2.Name = "button2";
            this.button2.Size = new System.Drawing.Size(100, 28);
            this.button2.TabIndex = 6;
            this.button2.Text = "Запуск";
            this.button2.UseVisualStyleBackColor = true;
            this.button2.Click += new System.EventHandler(this.button2_Click);
            // 
            // groupBox6
            // 
            this.groupBox6.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.groupBox6.Controls.Add(this.b_remove_player);
            this.groupBox6.Controls.Add(this.b_edit_player);
            this.groupBox6.Controls.Add(this.b_add_player);
            this.groupBox6.Controls.Add(this.lb_players);
            this.groupBox6.Location = new System.Drawing.Point(12, 48);
            this.groupBox6.Name = "groupBox6";
            this.groupBox6.Size = new System.Drawing.Size(451, 205);
            this.groupBox6.TabIndex = 8;
            this.groupBox6.TabStop = false;
            this.groupBox6.Text = "Список персонажей";
            // 
            // b_remove_player
            // 
            this.b_remove_player.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Right)));
            this.b_remove_player.Enabled = false;
            this.b_remove_player.Location = new System.Drawing.Point(343, 94);
            this.b_remove_player.Margin = new System.Windows.Forms.Padding(4);
            this.b_remove_player.Name = "b_remove_player";
            this.b_remove_player.Size = new System.Drawing.Size(100, 28);
            this.b_remove_player.TabIndex = 11;
            this.b_remove_player.Text = "Удалить";
            this.b_remove_player.UseVisualStyleBackColor = true;
            this.b_remove_player.Click += new System.EventHandler(this.b_remove_player_Click);
            // 
            // b_edit_player
            // 
            this.b_edit_player.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Right)));
            this.b_edit_player.Enabled = false;
            this.b_edit_player.Location = new System.Drawing.Point(344, 58);
            this.b_edit_player.Margin = new System.Windows.Forms.Padding(4);
            this.b_edit_player.Name = "b_edit_player";
            this.b_edit_player.Size = new System.Drawing.Size(100, 28);
            this.b_edit_player.TabIndex = 10;
            this.b_edit_player.Text = "Изменить";
            this.b_edit_player.UseVisualStyleBackColor = true;
            this.b_edit_player.Click += new System.EventHandler(this.b_edit_player_Click);
            // 
            // b_add_player
            // 
            this.b_add_player.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Right)));
            this.b_add_player.Location = new System.Drawing.Point(344, 22);
            this.b_add_player.Margin = new System.Windows.Forms.Padding(4);
            this.b_add_player.Name = "b_add_player";
            this.b_add_player.Size = new System.Drawing.Size(100, 28);
            this.b_add_player.TabIndex = 9;
            this.b_add_player.Text = "Добавить";
            this.b_add_player.UseVisualStyleBackColor = true;
            this.b_add_player.Click += new System.EventHandler(this.b_add_player_Click);
            // 
            // lb_players
            // 
            this.lb_players.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.lb_players.FormattingEnabled = true;
            this.lb_players.ItemHeight = 16;
            this.lb_players.Location = new System.Drawing.Point(7, 22);
            this.lb_players.Name = "lb_players";
            this.lb_players.Size = new System.Drawing.Size(328, 164);
            this.lb_players.TabIndex = 0;
            this.lb_players.SelectedIndexChanged += new System.EventHandler(this.lb_players_SelectedIndexChanged);
            // 
            // MainForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(8F, 16F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(479, 267);
            this.Controls.Add(this.groupBox6);
            this.Controls.Add(this.button2);
            this.Icon = ((System.Drawing.Icon)(resources.GetObject("$this.Icon")));
            this.KeyPreview = true;
            this.Margin = new System.Windows.Forms.Padding(4);
            this.MinimumSize = new System.Drawing.Size(497, 314);
            this.Name = "MainForm";
            this.Text = "TrainClient for Hunger Game Challenge";
            this.FormClosing += new System.Windows.Forms.FormClosingEventHandler(this.MainForm_FormClosing);
            this.Load += new System.EventHandler(this.MainForm_Load);
            this.KeyDown += new System.Windows.Forms.KeyEventHandler(this.MainForm_KeyDown);
            this.groupBox6.ResumeLayout(false);
            this.ResumeLayout(false);

        }

        #endregion
        private System.Windows.Forms.Button button2;
        private System.Windows.Forms.GroupBox groupBox6;
        private System.Windows.Forms.Button b_remove_player;
        private System.Windows.Forms.Button b_edit_player;
        private System.Windows.Forms.Button b_add_player;
        private System.Windows.Forms.ListBox lb_players;
    }
}

