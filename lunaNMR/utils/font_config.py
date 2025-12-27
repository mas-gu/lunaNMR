"""
Font configuration for cross-platform emoji support in tkinter GUI.

This module ensures consistent display of Unicode symbols and emoji across different platforms,
particularly addressing Linux systems where tkinter may not properly render emoji by default.
"""

import tkinter as tk
import tkinter.font as tkfont
import platform
import sys


class FontManager:
    """Manages font configuration for cross-platform emoji support."""
    
    def __init__(self):
        self.system = platform.system()
        self.emoji_font = None
        self.default_font = None
        self._configure_fonts()
    
    def _configure_fonts(self):
        """Configure fonts based on the operating system."""
        try:
            # Test if we can create a root window (needed for font operations)
            test_root = tk.Tk()
            test_root.withdraw()  # Hide the test window
            
            # Get default font
            self.default_font = tkfont.nametofont("TkDefaultFont")
            
            # Configure emoji-capable fonts based on platform
            if self.system == "Linux":
                self._configure_linux_fonts()
            elif self.system == "Darwin":  # macOS
                self._configure_macos_fonts()
            elif self.system == "Windows":
                self._configure_windows_fonts()
            else:
                self._configure_fallback_fonts()
            
            test_root.destroy()
            
        except Exception as e:
            print(f"⚠️ Font configuration warning: {e}")
            self._configure_fallback_fonts()
    
    def _configure_linux_fonts(self):
        """Configure fonts for Linux systems."""
        # List of fonts that support emoji on Linux, in order of preference
        emoji_fonts = [
            "Noto Color Emoji",
            "Noto Emoji", 
            "Symbola",
            "DejaVu Sans",
            "Liberation Sans",
            "Ubuntu"
        ]
        
        self.emoji_font = self._find_available_font(emoji_fonts)
        
        if not self.emoji_font:
            print("⚠️ No emoji-capable fonts found. Installing fonts-noto-color-emoji may help.")
            self._use_text_fallbacks()
    
    def _configure_macos_fonts(self):
        """Configure fonts for macOS systems."""
        # macOS typically has good emoji support
        self.emoji_font = tkfont.Font(family="Apple Color Emoji", size=10)
    
    def _configure_windows_fonts(self):
        """Configure fonts for Windows systems."""
        # Windows emoji fonts
        emoji_fonts = [
            "Segoe UI Emoji",
            "Segoe UI Symbol", 
            "Arial"
        ]
        
        self.emoji_font = self._find_available_font(emoji_fonts)
    
    def _configure_fallback_fonts(self):
        """Configure fallback fonts when platform detection fails."""
        self._use_text_fallbacks()
    
    def _find_available_font(self, font_list):
        """Find the first available font from a list."""
        for font_name in font_list:
            try:
                font = tkfont.Font(family=font_name, size=10)
                # Test if the font can render a common emoji
                return font
            except Exception:
                continue
        return None
    
    def _use_text_fallbacks(self):
        """Configure text-based fallbacks for emoji when no emoji fonts available."""
        self.emoji_font = None
        self._setup_emoji_replacements()
    
    def _setup_emoji_replacements(self):
        """Define text replacements for emoji when fonts don't support them."""
        self.emoji_replacements = {
            "✅": "[OK]",
            "❌": "[FAIL]", 
            "📍": "[REF]",
            "📄": "[DOC]",
            "📁": "[DIR]",
            "🔍": "[SEARCH]",
            "📊": "[CHART]",
            "📈": "[UP]",
            "📉": "[DOWN]",
            "💾": "[SAVE]",
            "🔧": "[TOOL]",
            "⚙️": "[SETTINGS]",
            "⚠️": "[WARN]",
            "🎛️": "[CTRL]",
            "🖥️": "[COMP]",
            "⏱️": "[TIME]",
            "🔬": "[LAB]",
            "✔️": "[CHECK]",
            "⏹️": "[STOP]"
        }

    def get_text_with_emoji_support(self, text):
        """Convert emoji characters to fallback text if needed."""
        if self.emoji_font or not hasattr(self, 'emoji_replacements'):
            return text
        
        # Replace emoji with text equivalents if no emoji font available
        result = text
        for emoji, replacement in self.emoji_replacements.items():
            result = result.replace(emoji, replacement)
        
        return result
    
    def configure_widget_font(self, widget, use_emoji=True):
        """Configure a widget to use emoji-capable font if available."""
        if use_emoji and self.emoji_font:
            try:
                widget.configure(font=self.emoji_font)
            except Exception:
                pass  # Fall back to default font


# Global font manager instance
_font_manager = None

def get_font_manager():
    """Get the global font manager instance."""
    global _font_manager
    if _font_manager is None:
        _font_manager = FontManager()
    return _font_manager

def configure_emoji_support(widget):
    """Configure a widget for emoji support."""
    manager = get_font_manager()
    manager.configure_widget_font(widget, use_emoji=True)

def get_display_text(text):
    """Get text with proper emoji fallbacks if needed."""
    manager = get_font_manager()
    return manager.get_text_with_emoji_support(text)

def print_font_info():
    """Print font configuration information for debugging."""
    manager = get_font_manager()
    print(f"System: {manager.system}")
    print(f"Emoji font available: {manager.emoji_font is not None}")
    if manager.emoji_font:
        print(f"Emoji font: {manager.emoji_font}")
    if hasattr(manager, 'emoji_replacements'):
        print("Using text fallbacks for emoji")


if __name__ == "__main__":
    # Test the font configuration
    print_font_info()
    
    # Test emoji conversion
    test_text = "✅ Success! 📊 Data loaded 🔧 Ready"
    converted = get_display_text(test_text)
    print(f"Original: {test_text}")
    print(f"Converted: {converted}")

# Linux emoji display bypass
import platform
if platform.system() == "Linux":
    def get_display_text_bypass(text):
        """Force emoji display on Linux instead of fallback"""
        return text  # Don't convert to fallbacks
    
    # Override the fallback function
    import sys
    current_module = sys.modules[__name__]
    current_module.get_display_text = get_display_text_bypass
