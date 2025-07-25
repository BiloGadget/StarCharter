import tkinter as tk
from tkinter import messagebox, simpledialog, colorchooser, filedialog
import math, json, os, copy

STAR_SAVE_DIR = "Chart Saves"
DEFAULT_SCALE = 2      # pixels per degree
MIN_SCALE = 1
MAX_SCALE = 2000
OVERLAP_PIXELS = 3     # base overlap threshold

# Default color coding by spectral type
TYPE_COLORS = {
    "O": "#95d0fc", "B": "#aabfff",
    "A": "#cad7ff", "F": "#f8f7ff",
    "G": "#fff4ea", "K": "#ffd2a1",
    "M": "#ffcc6f", "": "#ffffff"
}


class StarChartApp:
    def __init__(self, root):
        self.root = root
        root.title("Star Chart")
        root.resizable(False, False)

        # State
        self.stars = []
        self.mode = "Equirectangular"
        self.scale = DEFAULT_SCALE
        self.offset_x = 0
        self.offset_y = 0
        self.max_labels = 10
        self.mag_factor = 1.0           # global magnitude size multiplier
        self.undo_stack = []
        self.selected = None

        self._build_ui()
        self._load_or_create_save_dir()
        self._draw()

    def _build_ui(self):
        # Top controls
        top = tk.Frame(self.root)
        top.pack(fill=tk.X, padx=4, pady=4)
        tk.Button(top, text="Save",  command=self._save).pack(side=tk.LEFT)
        tk.Button(top, text="Load",  command=self._load).pack(side=tk.LEFT)
        tk.Button(top, text="Undo",  command=self._undo).pack(side=tk.LEFT)
        self.mode_btn = tk.Button(top, text="Mode", command=self._toggle_mode)
        self.mode_btn.pack(side=tk.LEFT, padx=4)
        self.color_btn = tk.Button(top, text="Color", command=self._edit_color_map)
        self.color_btn.pack(side=tk.LEFT)

        # Main area: canvas + editor
        main = tk.Frame(self.root)
        main.pack()
        self.canvas = tk.Canvas(main, width=800, height=600, bg="black")
        self.canvas.pack(side=tk.LEFT)
        # bind pan/zoom/click
        self.canvas.bind("<Button-1>", self._on_left_click)
        self.canvas.bind("<Button-3>", self._start_pan)
        self.canvas.bind("<B3-Motion>", self._do_pan)
        self.canvas.bind("<MouseWheel>", self._on_wheel)
        self.canvas.bind("<Button-4>", self._zoom_in)
        self.canvas.bind("<Button-5>", self._zoom_out)

        # Star editor
        editor = tk.LabelFrame(main, text="Star Editor", padx=6, pady=6)
        editor.pack(side=tk.LEFT, padx=8)
        self.fields = {}
        for lbl in ("Name", "RA (hhmmss.ss)", "Dec (±ddmmss.ss)", "Mag", "Type"):
            tk.Label(editor, text=lbl).pack(anchor=tk.W)
            ent = tk.Entry(editor, width=20)
            ent.pack(anchor=tk.W, pady=2)
            self.fields[lbl] = ent
        tk.Label(editor, text="Description").pack(anchor=tk.W)
        self.desc_text = tk.Text(editor, width=20, height=4)
        self.desc_text.pack()
        btn_frame = tk.Frame(editor)
        btn_frame.pack(pady=6)
        tk.Button(btn_frame, text="Add Star",    command=self._add_star).pack(side=tk.LEFT)
        tk.Button(btn_frame, text="Update Star", command=self._update_star).pack(side=tk.LEFT, padx=4)

        # Bottom-left controls (zoom, labels, mag-size)
        bottom_left = tk.Frame(self.root)
        bottom_left.pack(side=tk.BOTTOM, anchor=tk.W, padx=10, pady=4)

        tk.Button(bottom_left, text="-", command=self._zoom_out).pack(side=tk.LEFT)
        tk.Button(bottom_left, text="+", command=self._zoom_in).pack(side=tk.LEFT, padx=(4,10))

        tk.Label(bottom_left, text="Max Labels:").pack(side=tk.LEFT)
        self.lab_entry = tk.Entry(bottom_left, width=4)
        self.lab_entry.insert(0, str(self.max_labels))
        self.lab_entry.pack(side=tk.LEFT)
        tk.Button(bottom_left, text="Update", command=self._update_max_labels)\
            .pack(side=tk.LEFT, padx=(4,10))

        tk.Label(bottom_left, text="Mag Size:").pack(side=tk.LEFT)
        self.mag_slider = tk.Scale(
            bottom_left,
            from_=0.01, to=5.0,
            resolution=0.01,
            orient=tk.HORIZONTAL,
            length=640,
            command=self._on_mag_scale
        )
        self.mag_slider.set(self.mag_factor)
        self.mag_slider.pack(side=tk.LEFT)

        # Bottom-right delete
        self.delete_btn = tk.Button(self.root, text="Delete", bg="lightgray",
                                    command=self._confirm_delete)
        self.delete_btn.place(x=950, y=655)

    # File Save/Load -----------------------------------------------------
    def _load_or_create_save_dir(self):
        if not os.path.isdir(STAR_SAVE_DIR):
            os.makedirs(STAR_SAVE_DIR)

    def _save(self):
        name = simpledialog.askstring("Save Chart", "Enter filename:")
        if not name: return
        path = os.path.join(STAR_SAVE_DIR, name + ".json")
        with open(path, "w") as f:
            json.dump(self.stars, f, indent=2)
        messagebox.showinfo("Saved", f"Chart saved to {path}")

    def _load(self):
        files = [f for f in os.listdir(STAR_SAVE_DIR) if f.endswith(".json")]
        if not files:
            messagebox.showwarning("No Saves", "No saved charts in directory.")
            return

        win = tk.Toplevel(self.root)
        win.title("Load Chart")
        for fn in files:
            def _load_file(f=fn):
                path = os.path.join(STAR_SAVE_DIR, f)
                with open(path) as src:
                    self._push_undo()
                    self.stars = json.load(src)
                win.destroy()
                self._draw()
            btn = tk.Button(win, text=fn, width=30, command=_load_file)
            btn.pack(pady=2)

    # Undo ---------------------------------------------------------------
    def _push_undo(self):
        state = {
            "stars": copy.deepcopy(self.stars),
            "max_labels": self.max_labels
        }
        self.undo_stack.append(state)
        if len(self.undo_stack) > 20:  # cap history
            self.undo_stack.pop(0)

    def _undo(self):
        if not self.undo_stack: return
        state = self.undo_stack.pop()
        self.stars = state["stars"]
        self.max_labels = state["max_labels"]
        self.lab_entry.delete(0, tk.END)
        self.lab_entry.insert(0, str(self.max_labels))
        self._draw()

    # Mode & Color -------------------------------------------------------
    def _toggle_mode(self):
        self.mode = "Circular" if self.mode == "Equirectangular" else "Equirectangular"
        self._draw()

    def _edit_color_map(self):
        win = tk.Toplevel(self.root)
        win.title("Edit Color Map")
        for t in TYPE_COLORS:
            def _choose_color(tt=t):
                c = colorchooser.askcolor()[1]
                if c:
                    TYPE_COLORS[tt] = c
                    self._draw()
                    self._edit_color_map()  # refresh
                    win.destroy()
            btn = tk.Button(win, text=f"{t or 'Blank'}: {TYPE_COLORS[t]}",
                            bg=TYPE_COLORS[t], command=lambda tt=t: _choose_color(tt))
            btn.pack(fill=tk.X, pady=2)

    # Star management ----------------------------------------------------
    def _parse_ra(self, s):
        # expects hhmmss.s
        try:
            h = float(s[0:2])
            m = float(s[2:4])
            sec = float(s[4:])
            deg = (h + m/60 + sec/3600) * 15
            return deg
        except:
            raise ValueError("Bad RA format")

    def _parse_dec(self, s):
        # expects ±ddmmss.s
        try:
            sign = 1 if s[0] == "+" else -1
            d = abs(float(s[1:3]))
            m = float(s[3:5])
            sec = float(s[5:])
            return sign * (d + m/60 + sec/3600)
        except:
            raise ValueError("Bad Dec format")

    def _add_star(self):
        try:
            ra = self._parse_ra(self.fields["RA (hhmmss.ss)"].get())
            dec = self._parse_dec(self.fields["Dec (±ddmmss.ss)"].get())
        except ValueError as e:
            messagebox.showerror("Parse Error", str(e))
            return
        mag = float(self.fields["Mag"].get() or 2)
        star = {
            "name": self.fields["Name"].get(),
            "ra": ra, "dec": dec, "mag": mag,
            "type": self.fields["Type"].get(),
            "desc": self.desc_text.get("1.0", tk.END).strip()
        }
        self._push_undo()
        self.stars.append(star)
        self._draw()

    def _update_star(self):
        if self.selected is None: return
        idx = self.selected
        try:
            ra = self._parse_ra(self.fields["RA (hhmmss.ss)"].get())
            dec = self._parse_dec(self.fields["Dec (±ddmmss.ss)"].get())
        except ValueError as e:
            messagebox.showerror("Parse Error", str(e))
            return
        self._push_undo()
        st = self.stars[idx]
        st.update({
            "name": self.fields["Name"].get(),
            "ra": ra, "dec": dec,
            "mag": float(self.fields["Mag"].get() or 2),
            "type": self.fields["Type"].get(),
            "desc": self.desc_text.get("1.0", tk.END).strip()
        })
        self._draw()

    def _confirm_delete(self):
        if self.selected is None: return
        if messagebox.askyesno("Delete Star", "Confirm deletion?"):
            self._push_undo()
            del self.stars[self.selected]
            self.selected = None
            self._draw()

    # Label limit update -------------------------------------------------
    def _update_max_labels(self):
        try:
            n = int(self.lab_entry.get())
            if n < 0: raise ValueError()
        except:
            messagebox.showerror("Input Error", "Enter a positive integer.")
            return
        self._push_undo()
        self.max_labels = n
        self._draw()

    # Pan & Zoom ---------------------------------------------------------
    def _start_pan(self, ev):
        self.pan_start = (ev.x, ev.y, self.offset_x, self.offset_y)

    def _do_pan(self, ev):
        sx, sy, ox, oy = self.pan_start
        dx = (ev.x - sx) / self.scale
        dy = (sy - ev.y) / self.scale
        self.offset_x = ox - dx
        self.offset_y = oy - dy
        self._draw()

    def _zoom_in(self, ev=None):
        self.scale = min(MAX_SCALE, self.scale * 1.2)
        self._draw()

    def _zoom_out(self, ev=None):
        self.scale = max(MIN_SCALE, self.scale / 1.2)
        self._draw()

    def _on_wheel(self, ev):
        if ev.delta > 0:
            self._zoom_in()
        else:
            self._zoom_out()

    # Click to select or add by click ------------------------------------
    def _on_left_click(self, ev):
        # Check for star hit
        hit = self._find_star_at(ev.x, ev.y)
        if hit is not None:
            self.selected = hit
            self._populate_fields(self.stars[hit])
            self.delete_btn.config(bg="red")
        else:
            # add star by click
            ra, dec = self._screen_to_coord(ev.x, ev.y)
            self.fields["RA (hhmmss.ss)"].delete(0, tk.END)
            self.fields["Dec (±ddmmss.ss)"].delete(0, tk.END)
            self.fields["RA (hhmmss.ss)"].insert(0, self._fmt_ra(ra))
            self.fields["Dec (±ddmmss.ss)"].insert(0, self._fmt_dec(dec))
            self.selected = None
            self.delete_btn.config(bg="lightgray")

    def _find_star_at(self, x, y):
        # find topmost star within radius
        for i, st in enumerate(self.stars):
            sx, sy = self._coord_to_screen(st["ra"], st["dec"])
            r = max(0.05, (0.01 - st["mag"])) * (self.scale / DEFAULT_SCALE)
            if (x - sx)**2 + (y - sy)**2 <= (r + 2)**2:
                return i
        return None

    def _populate_fields(self, st):
        self.fields["Name"].delete(0, tk.END)
        self.fields["Name"].insert(0, st["name"])
        self.fields["RA (hhmmss.ss)"].delete(0, tk.END)
        self.fields["RA (hhmmss.ss)"].insert(0, self._fmt_ra(st["ra"]))
        self.fields["Dec (±ddmmss.ss)"].delete(0, tk.END)
        self.fields["Dec (±ddmmss.ss)"].insert(0, self._fmt_dec(st["dec"]))
        self.fields["Mag"].delete(0, tk.END)
        self.fields["Mag"].insert(0, str(st["mag"]))
        self.fields["Type"].delete(0, tk.END)
        self.fields["Type"].insert(0, st["type"])
        self.desc_text.delete("1.0", tk.END)
        self.desc_text.insert(tk.END, st["desc"])

    # Coordinate transforms ----------------------------------------------
    def _coord_to_screen(self, ra_deg, dec_deg):
        cx, cy = 400, 300
        if self.mode == "Equirectangular":
            x = cx + (ra_deg - self.offset_x) * self.scale
            y = cy - (dec_deg - self.offset_y) * self.scale
        else:
            # Circular projection
            r = (90 - dec_deg) * self.scale
            ang = math.radians(ra_deg)
            x = cx + r * math.cos(ang)
            y = cy + r * math.sin(ang)
        return x, y

    def _screen_to_coord(self, x, y):
        cx, cy = 400, 300
        if self.mode == "Equirectangular":
            ra = (x - cx) / self.scale + self.offset_x
            dec = (cy - y) / self.scale + self.offset_y
        else:
            dx, dy = x - cx, y - cy
            r = math.hypot(dx, dy)
            dec = 90 - r / self.scale
            ang = math.atan2(dy, dx)
            ra = (math.degrees(ang)) % 360
        return ra, dec

    def _fmt_ra(self, deg):
        h = int(deg / 15)
        m = int((deg / 15 - h) * 60)
        s = (deg / 15 - h - m/60) * 3600
        return f"{h:02d}{m:02d}{s:05.2f}"

    def _fmt_dec(self, deg):
        sign = "+" if deg >= 0 else "-"
        d = int(abs(deg))
        m = int((abs(deg) - d) * 60)
        s = (abs(deg) - d - m/60) * 3600
        return f"{sign}{d:02d}{m:02d}{s:05.2f}"

    # Drawing ----------------------------------------------------------------
    def _draw(self):
        self.canvas.delete("all")
        self._draw_grid()
        self._draw_stars()
        self.delete_btn.config(bg="lightgray")

    def _draw_grid(self):
        w, h = 800, 600
        cx, cy = w//2, h//2
        if self.mode == "Equirectangular":
            # vertical RA lines every 15° (1h)
            for ra in range(0, 361, 15):
                x, _ = self._coord_to_screen(ra, 0)
                self.canvas.create_line(x, 0, x, h, fill="#333333")
                self.canvas.create_text(x+2, 10, text=f"{ra//15}h", fill="#888")

            # horizontal Dec every 15°
            for dec in range(-90, 91, 15):
                _, y = self._coord_to_screen(0, dec)
                self.canvas.create_line(0, y, w, y, fill="#333333")
                self.canvas.create_text(10, y-2, text=f"{dec}°", fill="#888", anchor=tk.W)
        else:
            # circular border
            rmax = min(cx, cy) - 5
            self.canvas.create_oval(cx-rmax, cy-rmax, cx+rmax, cy+rmax, outline="#333333")
            # radial hour lines
            for h in range(24):
                ang = math.radians(h * 15)
                x1 = cx + rmax * math.cos(ang)
                y1 = cy + rmax * math.sin(ang)
                x2 = cx + (rmax - 10) * math.cos(ang)
                y2 = cy + (rmax - 10) * math.sin(ang)
                self.canvas.create_line(x1, y1, x2, y2, fill="#333333")
                lx = cx + (rmax+12) * math.cos(ang)
                ly = cy + (rmax+12) * math.sin(ang)
                self.canvas.create_text(lx, ly, text=f"{h}h", fill="#888")

    def _on_mag_scale(self, val):
        self.mag_factor = float(val)
        self._draw()

    # re-define _draw_stars to apply mag_factor
    def _draw_stars(self):
        onscreen = []
        for st in self.stars:
            x, y = self._coord_to_screen(st["ra"], st["dec"])
            if 0 <= x <= 800 and 0 <= y <= 600:
                onscreen.append((st["mag"], st))
        onscreen.sort(key=lambda x: x[0])

        displayed = []
        for mag, st in onscreen:
            x, y = self._coord_to_screen(st["ra"], st["dec"])
            base_r = max(2, (10 - mag))
            r = base_r * self.mag_factor * (self.scale / DEFAULT_SCALE)
            overlap = False
            for xx, yy, rr in displayed:
                if (x - xx)**2 + (y - yy)**2 < (r + rr + OVERLAP_PIXELS)**2:
                    overlap = True
                    break
            if overlap: continue
            displayed.append((x, y, r))
            color = TYPE_COLORS.get(st["type"], "#fff")
            self.canvas.create_oval(x-r, y-r, x+r, y+r, fill=color, outline="")

        for mag, st in onscreen[:self.max_labels]:
            x, y = self._coord_to_screen(st["ra"], st["dec"])
            self.canvas.create_text(x+5, y-5, text=st["name"],
                                    fill="#fff", anchor=tk.SW)

# Main ------------------------------
if __name__ == "__main__":
    root = tk.Tk()
    app = StarChartApp(root)
    root.mainloop()
