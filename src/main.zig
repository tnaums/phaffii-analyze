const std = @import("std");
const fasta = @import("fasta");
const Io = std.Io;

const promoter = enum { aox1, gap, unknown };
const secretion = enum { alpha, ost, cytoplasmic };

const Komagataella = struct {
    plasmid: fasta.DNA,
    promoter: promoter,
    coding: fasta.DNA,
    secretion: secretion,
    recombinant: fasta.Protein,

    pub fn init(io: Io, allocator: std.mem.Allocator, filepath: []const u8) !Komagataella {
        var dna = try parseSingleDNA(io, allocator, filepath);
        try dna.addTranslation(allocator);
        const p = checkPromoter(dna);
        var codingRegion = try getCodingRegion(allocator, dna, p);
        try codingRegion.addTranslation(allocator);
        const s = determineSecretion(codingRegion);
        const r = try recombinantProtein(allocator, codingRegion, s);
        return Komagataella{
            .plasmid = dna,
            .promoter = p,
            .coding = codingRegion,
            .secretion = s,
            .recombinant = r,
        };
    }

    pub fn deinit(self: Komagataella, allocator: std.mem.Allocator) void {
        self.plasmid.deinit(allocator);
        self.coding.deinit(allocator);
        self.recombinant.deinit(allocator);
    }

    fn checkPromoter(dna: fasta.DNA) promoter {
        const aox1Start = std.mem.indexOf(u8, dna.sequence, "AGATCTAACATCCAAA");
        const aox1End = std.mem.indexOf(u8, dna.sequence, "ATTCGAAACGA") orelse 0;

        const gapStart = std.mem.indexOf(u8, dna.sequence, "AGATCTTTTTTGTAGAAATG");
        const gapEnd = std.mem.indexOf(u8, dna.sequence, "TATTTCAATCAATTGAACAAC") orelse 0;

        if (aox1Start) |value| {
            if (value == 0 and aox1End == 930) {
                return promoter.aox1;
            }
        }
        if (gapStart) |value| {
            if (value == 0 and gapEnd == 459) {
                return promoter.gap;
            }
        }
        return promoter.unknown;
    }

    fn getCodingRegion(allocator: std.mem.Allocator, dna: fasta.DNA, p: promoter) !fasta.DNA {
        switch(p) {
            .aox1 => {
                const codingEnd = std.mem.indexOf(u8, dna.sequence, "GTTTGTAGCCTTAGA") orelse dna.sequence.len;
                return try fasta.DNA.init(allocator, "CodingRegion", dna.sequence[940..codingEnd]);
            },
            .gap => {
                const codingEnd = std.mem.indexOf(u8, dna.sequence, "GTTTTAGCCTTAGAC") orelse dna.sequence.len;                
                return try fasta.DNA.init(allocator, "CodingRegion", dna.sequence[492..codingEnd]);
            },
            .unknown => {unreachable;},
        }
    }

    fn determineSecretion(dna: fasta.DNA) secretion {
        if  (std.mem.startsWith(u8, dna.sequence, "ATGAGATTTCCTTCA")) {
            return secretion.alpha;
        }
        if (std.mem.startsWith(u8, dna.sequence, "ATGAGGCAGGTTTGG")) {
            return secretion.ost;
        }
        return secretion.cytoplasmic;
    }

    fn recombinantProtein(allocator: std.mem.Allocator, codingRegion: fasta.DNA, s: secretion) !fasta.Protein {
        var fullP: []u8 = undefined;
        std.debug.print("secretion value: {any}\n", .{s});
        if (codingRegion.translation) |value| {
            fullP = value[0];
        }
        var parts = std.mem.tokenizeScalar(u8, fullP, '*');
        var shortLeft = parts.next() orelse fullP;

        std.debug.print("{s}\n", .{fullP});
        std.debug.print("{s}\n", .{shortLeft});

        switch (s) {
            .alpha => {
                return try fasta.Protein.init(allocator, "mature", shortLeft[89..]);
            },

            .ost => {
                return try fasta.Protein.init(allocator, "mature", shortLeft[92..]);
            },
            .cytoplasmic => {
                return try fasta.Protein.init(allocator, "mature", shortLeft);
            },
        }
    }
};

fn isSecreted(k: Komagataella) bool {
    if (std.mem.eql(u8, k.plasmid.sequence[940..949], "ATGAGATTT") and
        (std.mem.eql(u8, k.plasmid.sequence[1198..1207], "GCTGAAGCT")))
    {
        return true;
    }
    return false;
}

pub fn parseSingleDNA(io: Io, allocator: std.mem.Allocator, filepath: []const u8) !fasta.DNA {
    // open file
    const file = try std.Io.Dir.cwd().openFile(io, filepath, .{});
    defer file.close(io);

    const state = enum { inHeader, inSequence };
    var myState: ?state = null;
    // initialize ArrayLists
    var header = std.ArrayList(u8).empty;
    defer header.deinit(allocator);
    var sequence = std.ArrayList(u8).empty;
    defer sequence.deinit(allocator);

    var buf: [1024]u8 = undefined;
    // parse out the info
    while (true) {
        const n = file.readStreaming(io, &.{&buf}) catch |err| {
            if (err == error.EndOfStream) break;
            return err;
        };
        if (n == 0) {
            if (header.items.len == 0) break;
            break;
        }
        var i: u16 = 0;
        while (i < n) : (i += 1) {
            if (myState) |s| {
                switch (s) {
                    .inHeader => {
                        if (buf[i] == '\n') {
                            myState = state.inSequence;
                            continue;
                        }
                        try header.append(allocator, buf[i]);
                    },
                    .inSequence => {
                        if (buf[i] == '>') { // if a second sequence begins, just return the first one
                            return try fasta.DNA.init(allocator, header.items, sequence.items);
                        }
                        if (buf[i] != '\n') {
                            try sequence.append(allocator, std.ascii.toUpper(buf[i]));
                        }
                    },
                }
            } else {
                if (buf[i] == '>') {
                    myState = state.inHeader;
                    continue;
                }
            }
        }
    }
    return try fasta.DNA.init(allocator, header.items, sequence.items);
}

pub fn main(init: std.process.Init) !void {
    const stdout = std.Io.File.stdout();

    const args = try init.minimal.args.toSlice(init.arena.allocator());
    if (args.len != 2) {
        std.debug.print("Usage: {s} <filename.fa>\n", .{args[0]});
        return;
    }

    const filepath = args[1];

    var k = try Komagataella.init(init.io, init.gpa, filepath);
    defer k.deinit(init.gpa);

    const fs = try std.fmt.allocPrint(init.gpa, "{f}\n", .{k.plasmid});
    defer init.gpa.free(fs);
    try stdout.writeStreamingAll(init.io, fs);

//    try k.plasmid.mapDNA(init.gpa, init.io, stdout);
    switch (k.promoter) {
        .aox1 => {
            try stdout.writeStreamingAll(init.io, "Promoter is inducible aox1.\n");
        },
        .gap => {
            try stdout.writeStreamingAll(init.io, "Promoter is constitutive gap.\n");
        },
        .unknown => {
            try stdout.writeStreamingAll(init.io, "Unknown promoter type.\n");
        },
    }
    switch (k.secretion) {
        .alpha => {
            try stdout.writeStreamingAll(init.io, "Secreted protein with alpha factor.\n");
        },
        .ost => {
            try stdout.writeStreamingAll(init.io, "Secreted protein with ost1.\n");
        },
        .cytoplasmic => {
            try stdout.writeStreamingAll(init.io, "Possible cytoplasmic protein.\n");
        }
    }
    const mp = try std.fmt.allocPrint(init.gpa, "{f}\n", .{k.recombinant});
    defer init.gpa.free(mp);
    try stdout.writeStreamingAll(init.io, mp);

}
