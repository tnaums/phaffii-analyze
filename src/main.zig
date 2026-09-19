const std = @import("std");
const fasta = @import("fasta");
const Io = std.Io;

const Komagataella = struct {
    plasmid: fasta.DNA,
};

pub fn main(init: std.process.Init) !void {
    const stdout = std.Io.File.stdout();
    
    const args = try init.minimal.args.toSlice(init.arena.allocator());
    if (args.len != 2) {
        std.debug.print("Usage: {s} <filename.fa>\n", .{args[0]});
        return;
    }

    const filepath = args[1];
    const file = try std.Io.Dir.cwd().openFile(init.io, filepath, .{});
    defer file.close(init.io);

    // Starting DNA queue
    var queue: std.Io.Queue(fasta.DNA) = .init(&.{});
    var producer_task = try init.io.concurrent(fasta.parseDNA, .{ init.io, init.gpa, &queue, file });
    defer producer_task.cancel(init.io) catch {};

    var myPlasmid = try queue.getOne(init.io);
    defer myPlasmid.deinit(init.gpa);

    const k = Komagataella{
        .plasmid = myPlasmid,
    };

    const fs = try std.fmt.allocPrint(init.gpa, "{f}\n", .{k.plasmid});
    defer init.gpa.free(fs);
    try stdout.writeStreamingAll(init.io, fs);

    try myPlasmid.mapDNA(init.gpa, init.io, stdout);
    
}
